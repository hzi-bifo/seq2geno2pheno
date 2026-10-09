#!/usr/bin/env nextflow
/*
 * Seq2Geno2Pheno v3 -- assemblies in, evaluated phenotype models (or a trained model) out.
 *
 *   nextflow run hzi-bifo/seq2geno2pheno -profile singularity \
 *       --input samples.tsv --preset noleak --outdir results
 *
 * A preset (src/s2g2p/presets/<name>.yml) fixes every method choice; this file only wires the
 * steps. `features.fit_scope` decides the shape of the run:
 *
 *   cohort  one Roary pangenome over the whole sample sheet, compressed once (hu2024, g2p)
 *   fold    one pangenome per (phenotype, outer fold) from its training genomes, the held-out
 *           genomes placed into it (noleak, fast)
 *
 * --mode train fits a final model per phenotype on all labelled genomes (fold-scope presets
 * only); --mode predict places new genomes into that model's pangenome and predicts them.
 */
nextflow.enable.dsl = 2

include { PROKKA; ROARY; KMERS; GPA; INDELS; COMPRESS; SPLIT; PLACE; EVALUATE; TRAIN; PREDICT } from './modules/local/s2g2p'

def loadPreset(name) {
    def f = file("${projectDir}/src/s2g2p/presets/${name}.yml")
    if (!f.exists()) {
        f = file(name)
    }
    if (!f.exists()) {
        error "unknown preset '${name}' (see: s2g2p presets)"
    }
    return new org.yaml.snakeyaml.Yaml().load(f.text)
}

workflow {
    if (!params.input) {
        error "--input <sample sheet> is required"
    }
    def preset = loadPreset(params.preset)
    def feats = preset.features
    log.info "Seq2Geno2Pheno v3 | preset ${preset.name} (${preset.status}) | mode ${params.mode} | fit_scope ${feats.fit_scope}"

    sheet = file(params.input, checkIfExists: true)
    folds = params.folds ? file(params.folds.replace('{phenotype}', '*')) : []

    samples = Channel.fromPath(sheet)
        .splitCsv(header: true, sep: sheet.name.endsWith('.csv') ? ',' : '\t')
        // relative assembly paths resolve against the sheet's directory, as in s2g2p.io
        .map { row -> tuple(row.sample, file(row.assembly.startsWith('/') ? row.assembly : "${sheet.parent}/${row.assembly}", checkIfExists: true)) }

    annotations = PROKKA(samples)                       // [sample, dir <sample>/]
    annot_dirs = annotations.map { it[1] }.collect()
    kmers = KMERS(sheet)

    if (params.mode == 'evaluate' && feats.fit_scope == 'cohort') {
        all_ids = samples.map { it[0] }.collectFile(name: 'cohort.txt', newLine: true)
        pg = ROARY(all_ids.map { tuple('cohort', it) }, annot_dirs, feats.roary_args)
        gpa = GPA(pg)
        ind = INDELS(gpa.join(pg), annot_dirs, feats.indel_min_isolates, feats.indel_legacy_rules)
        tables = COMPRESS(gpa.join(ind.tables))          // [label, gpa_NONRDNT.tsv, indel_NONRDNT.tsv]
        EVALUATE(sheet, folds, kmers, tables.map { it[1..2] }.collect(), params.preset, 'cohort')
    }
    else if (params.mode == 'evaluate' && feats.fit_scope == 'fold') {
        splits = SPLIT(sheet, folds)
        jobs = splits.lists.flatten()
            .filter { it.name.endsWith('.train.txt') }
            .map { f ->
                def ph = f.parent.name
                def k = (f.name =~ /fold_(\d+)\.train\.txt/)[0][1]
                tuple("${ph}__fold_${k}".toString(), f, file(f.toString().replace('.train.txt', '.test.txt')))
            }
        pg = ROARY(jobs.map { tuple(it[0], it[1]) }, annot_dirs, feats.roary_args)
        gpa = GPA(pg)
        ind = INDELS(gpa.join(pg), annot_dirs, feats.indel_min_isolates, feats.indel_legacy_rules)
        placed = PLACE(gpa.join(ind.tables).join(ind.aln).join(pg).join(jobs.map { tuple(it[0], it[2]) }), annot_dirs)
        EVALUATE(sheet, folds, kmers, placed.map { it[1] }.collect(), params.preset, 'fold')
    }
    else if (params.mode == 'train') {
        if (!preset.predict) {
            error "preset ${preset.name} cannot make a model for unseen genomes (predict: false); use noleak"
        }
        all_ids = samples.map { it[0] }.collectFile(name: 'cohort.txt', newLine: true)
        pg = ROARY(all_ids.map { tuple('final', it) }, annot_dirs, feats.roary_args)
        gpa = GPA(pg)
        ind = INDELS(gpa.join(pg), annot_dirs, feats.indel_min_isolates, feats.indel_legacy_rules)
        TRAIN(sheet, folds, kmers, gpa.join(ind.tables).join(ind.aln).join(pg), params.preset)
    }
    else if (params.mode == 'predict') {
        if (!params.model_dir) {
            error "--mode predict needs --model_dir (the outdir of a --mode train run)"
        }
        PREDICT(samples.map { it[0] }.collectFile(name: 'query.txt', newLine: true), annot_dirs, kmers,
                file(params.model_dir, checkIfExists: true))
    }
    else {
        error "--mode must be evaluate, train or predict"
    }
}
