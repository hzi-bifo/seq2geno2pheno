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
 * --mode folds writes outer folds per phenotype (--fold_method random | phylogeny | kma) for a
 * later --folds; --mode install creates the tool environments for -profile conda (envs/install.sh).
 */
nextflow.enable.dsl = 2

include { INSTALL; PROKKA; ROARY; KMERS; GPA; INDELS; COMPRESS; SPLIT; NJ_TREE; KMA_CLUSTER; FOLDS; PLACE; EVALUATE; TRAIN; PREDICT } from './modules/local/s2g2p'

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

def phenotypeColumns(sheet) {
    def sep = sheet.name.endsWith('.csv') ? ',' : '\t'
    return sheet.readLines()[0].split(sep).collect { c -> c.trim() }.findAll { c -> !(c in ['sample', 'assembly', 'group']) }
}

workflow {
    // Nextflow's strict syntax (the default from 26.x) allows no `return` in a workflow body, so
    // install is the first branch of one if/else chain rather than an early exit.
    if (params.mode == 'install') {
        INSTALL().view()
    }
    else {
        run_pipeline()
    }
}

workflow run_pipeline {
    main:
    if (workflow.profile.tokenize(',').contains('conda') && !workflow.stubRun
            && !file("${params.env_dir}/.installed").exists()) {
        error "-profile conda: no environments under --env_dir ${params.env_dir}; install them first " +
              "with `bash ${projectDir}/envs/install.sh --prefix ${params.env_dir}` or --mode install"
    }
    if (!params.input) {
        error "--input <sample sheet> is required"
    }
    if (workflow.workDir == workflow.launchDir.resolve('work') && !workflow.stubRun) {
        log.warn "No -work-dir given: Prokka's annotations, every Roary pangenome and the alignments live only in " +
                 "${workflow.workDir}, on whatever storage that is. Choose it with -work-dir, and keep it until the run is " +
                 "checked: -resume needs it."
    }
    def preset = loadPreset(params.preset)
    def feats = preset.features
    log.info "Seq2Geno2Pheno v3 | preset ${preset.name} (${preset.status}) | mode ${params.mode} | fit_scope ${feats.fit_scope}"

    def sheet = file(params.input, checkIfExists: true)
    def folds = params.folds ? file(params.folds.replace('{phenotype}', '*')) : []

    def samples = channel.fromPath(sheet)
        .splitCsv(header: true, sep: sheet.name.endsWith('.csv') ? ',' : '\t')
        // relative assembly paths resolve against the sheet's directory, as in s2g2p.io
        .map { row -> tuple(row.sample, file(row.assembly.startsWith('/') ? row.assembly : "${sheet.parent}/${row.assembly}", checkIfExists: true)) }

    if (params.mode == 'folds') {
        // Folds only: nothing below (k-mers, features, models) runs. FOLDS writes
        // <outdir>/folds/<phenotype>_<method>_cv.json, the input of a later --folds.
        def method = params.fold_method
        if (!(method in ['random', 'phylogeny', 'kma'])) {
            error "--fold_method must be random, phylogeny or kma"
        }
        if (method != 'random' && !workflow.stubRun && !workflow.profile.tokenize(',').contains('conda')) {
            error "--fold_method ${method} needs -profile conda: no container image carries R/phangorn or kma_clustering yet"
        }
        def tool = method == 'phylogeny' ? 'nj_tree' : 'kma_clustering'
        if (method != 'random' && !workflow.stubRun && !(method == 'phylogeny' && params.tree)
                && !file("${params.env_dir}/bin/${tool}").exists()) {
            error "--fold_method ${method}: ${params.env_dir} has no ${tool} (an install from before --mode folds); " +
                  "rerun `bash ${projectDir}/envs/install.sh --prefix ${params.env_dir}`, which adds only the folds environment"
        }
        if (method == 'phylogeny') {
            def tree
            if (params.tree) {
                tree = channel.value(file(params.tree, checkIfExists: true))
            }
            else {
                if (!feats.roary_args.tokenize(' ').contains('-e')) {
                    error "--fold_method phylogeny needs Roary's core-gene alignment, which preset ${preset.name} skips (-e --mafft); use noleak or hu2024, or pass --tree"
                }
                // The same cohort pangenome as hu2024's (label, genome list and arguments), so a
                // later hu2024 run with -resume reuses it.
                def fold_annot = PROKKA(samples).map { s -> s[1] }.collect()
                def cohort_ids = samples.map { s -> s[0] }.collectFile(name: 'cohort.txt', newLine: true)
                tree = NJ_TREE(ROARY(cohort_ids.map { f -> tuple('cohort', f) }, fold_annot, feats.roary_args))
            }
            FOLDS(sheet, tree, [], method)
        }
        else if (method == 'kma') {
            def phenos = params.phenotype == 'all' ? phenotypeColumns(sheet) : [params.phenotype]
            def asm = samples.toList().map { l -> tuple(l.collect { s -> s[0] }, l.collect { s -> s[1] }) }
            FOLDS(sheet, [], KMA_CLUSTER(sheet, asm, phenos).collect(), method)
        }
        else {
            FOLDS(sheet, [], [], method)
        }
    }
    else {
        def annotations = PROKKA(samples)                   // [sample, dir <sample>/]
        def annot_dirs = annotations.map { s -> s[1] }.collect()
        def kmers = KMERS(samples.toList().map { l -> tuple(l.collect { s -> s[0] }, l.collect { s -> s[1] }) })

        if (params.mode == 'evaluate' && feats.fit_scope == 'cohort') {
            def all_ids = samples.map { s -> s[0] }.collectFile(name: 'cohort.txt', newLine: true)
            def pg = ROARY(all_ids.map { f -> tuple('cohort', f) }, annot_dirs, feats.roary_args)
            def gpa = GPA(pg)
            def ind = INDELS(gpa.join(pg), annot_dirs, feats.indel_min_isolates, feats.indel_legacy_rules)
            def tables = COMPRESS(gpa.join(ind.tables))      // [label, gpa_NONRDNT.tsv, indel_NONRDNT.tsv]
            EVALUATE(sheet, folds, kmers, tables.map { t -> t[1..2] }.collect(), params.preset, 'cohort', preset.ml_env ?: 'current')
        }
        else if (params.mode == 'evaluate' && feats.fit_scope == 'fold') {
            def splits = SPLIT(sheet, folds)
            def jobs = splits.lists.flatten()
                .filter { f -> f.name.endsWith('.train.txt') }
                .map { f ->
                    def ph = f.parent.name
                    def k = (f.name =~ /fold_(\d+)\.train\.txt/)[0][1]
                    tuple("${ph}__fold_${k}".toString(), f, file(f.toString().replace('.train.txt', '.test.txt')))
                }
            def pg = ROARY(jobs.map { j -> tuple(j[0], j[1]) }, annot_dirs, feats.roary_args)
            def gpa = GPA(pg)
            def ind = INDELS(gpa.join(pg), annot_dirs, feats.indel_min_isolates, feats.indel_legacy_rules)
            def placed = PLACE(gpa.join(ind.tables).join(ind.aln).join(pg).join(jobs.map { j -> tuple(j[0], j[2]) }), annot_dirs)
            EVALUATE(sheet, folds, kmers, placed.map { t -> t[1] }.collect(), params.preset, 'fold', preset.ml_env ?: 'current')
        }
        else if (params.mode == 'train') {
            if (!preset.predict) {
                error "preset ${preset.name} cannot make a model for unseen genomes (predict: false); use noleak"
            }
            def all_ids = samples.map { s -> s[0] }.collectFile(name: 'cohort.txt', newLine: true)
            def pg = ROARY(all_ids.map { f -> tuple('final', f) }, annot_dirs, feats.roary_args)
            def gpa = GPA(pg)
            def ind = INDELS(gpa.join(pg), annot_dirs, feats.indel_min_isolates, feats.indel_legacy_rules)
            TRAIN(sheet, folds, kmers, gpa.join(ind.tables).join(ind.aln).join(pg), annot_dirs, params.preset)
        }
        else if (params.mode == 'predict') {
            if (!params.model_dir) {
                error "--mode predict needs --model_dir (the outdir of a --mode train run)"
            }
            PREDICT(samples.map { s -> s[0] }.collectFile(name: 'query.txt', newLine: true), annot_dirs, kmers,
                    file(params.model_dir, checkIfExists: true))
        }
        else {
            error "--mode must be evaluate, train, predict, folds or install"
        }
    }
}
