/*
 * Seq2Geno2Pheno v3 processes. Each external tool runs from its own pinned environment: the
 * wrappers envs/install.sh writes (-profile conda), or an image (params.containers, container
 * profiles). Every Python step is the `s2g2p` CLI of this repository.
 *
 * Annotation directories are staged under their sample names, so `--annotations .` finds
 * <sample>/<sample>.{gff,faa,ffn} in the task directory.
 */

process INSTALL {
    // The environments of -profile conda. Runs where the executor puts it: that node needs
    // internet access, and --env_dir must be on storage the compute nodes share.
    label 'process_medium'

    output:
    stdout

    script:
    """
    bash ${projectDir}/envs/install.sh --prefix ${params.env_dir}
    """

    stub:
    """
    echo install stub
    """
}

process PROKKA {
    // Seq2Geno's create_gff rule, flag for flag (Prokka 1.14.6 in the published env).
    tag "${sample}"
    label 'process_medium'
    container params.containers.prokka
    publishDir "${params.outdir}/annotations", mode: params.publish_dir_mode, enabled: params.publish_annotations

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}", type: 'dir')

    script:
    """
    prokka --locustag ${sample} --prefix ${sample} --force --cpus ${task.cpus} \\
        --metagenome --compliant --outdir ${sample} ${fasta}
    """

    stub:
    """
    mkdir ${sample}
    touch ${sample}/${sample}.gff ${sample}/${sample}.faa ${sample}/${sample}.ffn
    """
}

process KMERS {
    // Assemblies are staged (asm/<n>/<file>), never read from the sheet's host paths: those
    // are not mounted inside the container.
    label 'process_medium'
    container params.containers.s2g2p
    publishDir "${params.outdir}/features", mode: params.publish_dir_mode

    input:
    tuple val(ids), path(fastas, stageAs: 'asm/?/*')

    output:
    path 'kmers.tsv'

    script:
    def rows = [ids, fastas instanceof List ? fastas : [fastas]].transpose().collect { s, f -> "${s}\t${f}" }.join('\n')
    """
    printf 'sample\tassembly\n${rows}\n' > assemblies.tsv
    s2g2p kmers --samplesheet assemblies.tsv --k ${params.kmer_k} --threads ${task.cpus} --out kmers.tsv
    """

    stub:
    """
    touch kmers.tsv
    """
}

process ROARY {
    // Seq2Geno's roary rule: Roary 3.8.2 (Seq2Geno's bundled build, its own BLAST 2.4.0) on the
    // Prokka GFFs of the genomes listed in `ids`.
    tag "${label}"
    label 'process_high'
    container params.containers.roary

    input:
    tuple val(label), path(ids)
    path annotations
    val roary_args

    output:
    tuple val(label), path('roary', type: 'dir')

    script:
    """
    mkdir gff
    while read -r s; do
        [ -n "\$s" ] && ln -s "\$(readlink -f \$s/\$s.gff)" gff/\$s.gff
    done < ${ids}
    n=\$(ls gff | wc -l)
    [ "\$n" -ge 2 ] || { echo "need at least two genomes, got \$n" >&2; exit 1; }
    roary -f roary -v gff/*.gff -p ${task.cpus} ${roary_args}
    # Roary exits 0 when a step it shells out to fails (e.g. a program missing from PATH), and
    # just leaves no tables; stop here rather than one step later with a less useful error.
    for f in gene_presence_absence.csv gene_presence_absence.Rtab clustered_proteins; do
        [ -s roary/\$f ] || { echo "Roary wrote no roary/\$f; its log is above" >&2; exit 1; }
    done
    """

    stub:
    """
    mkdir roary
    touch roary/gene_presence_absence.csv roary/clustered_proteins
    """
}

process GPA {
    tag "${label}"
    label 'process_low'
    container params.containers.s2g2p

    input:
    tuple val(label), path(roary)

    output:
    tuple val(label), path('gpa.tsv')

    script:
    """
    s2g2p gpa --roary-csv ${roary}/gene_presence_absence.csv --out gpa.tsv
    """

    stub:
    """
    touch gpa.tsv
    """
}

process INDELS {
    // Seq2Geno's indel_* rules: mafft per family, msa2vcf, vcf2indel (ported), the table.
    tag "${label}"
    label 'process_high'
    container params.containers.tools

    input:
    tuple val(label), path(gpa), path(roary)
    path annotations
    val min_isolates
    val legacy_rules

    output:
    tuple val(label), path('indel/indel.tsv'), emit: tables
    tuple val(label), path('indel/aln', type: 'dir'), emit: aln

    script:
    def legacy = legacy_rules ? '' : '--no-legacy-rules'
    """
    s2g2p indels --roary-dir ${roary} --annotations . --gpa ${gpa} --out-dir indel \\
        --min-isolates ${min_isolates} ${legacy} --threads ${task.cpus}
    """

    stub:
    """
    mkdir -p indel/aln
    touch indel/indel.tsv
    """
}

process COMPRESS {
    // featCompress over the whole cohort (cohort-scope presets only).
    tag "${label}"
    label 'process_low'
    container params.containers.s2g2p
    publishDir "${params.outdir}/features", mode: params.publish_dir_mode

    input:
    tuple val(label), path(gpa), path(indel)

    output:
    tuple val(label), path('gpa_NONRDNT.tsv'), path('indel_NONRDNT.tsv')

    script:
    """
    s2g2p compress --input ${gpa} --out gpa_NONRDNT.tsv --groups gpa_GROUPS.txt
    s2g2p compress --input ${indel} --out indel_NONRDNT.tsv --groups indel_GROUPS.txt
    """

    stub:
    """
    touch gpa_NONRDNT.tsv indel_NONRDNT.tsv
    """
}

process SPLIT {
    label 'process_single'
    container params.containers.s2g2p
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path sheet
    path folds

    output:
    path 'splits/*/fold_*.txt', emit: lists
    path 'splits/manifest.tsv'

    script:
    def f = params.folds ? "--folds '${file(params.folds).name}'" : "--n-folds ${params.n_folds} --seed ${params.seed}"
    """
    s2g2p split --samplesheet ${sheet} --phenotype ${params.phenotype} ${f} --out-dir splits
    """

    stub:
    """
    mkdir -p splits/drug
    for k in 0 1; do touch splits/drug/fold_\$k.train.txt splits/drug/fold_\$k.test.txt; done
    touch splits/manifest.tsv
    """
}

process NJ_TREE {
    // Phylogeny-aware folds: phangorn NJ on Hamming distances of Roary's core-gene alignment
    // (src/s2g2p/R/nj_tree.R; AMR_benchmarking src/cv_folds/phylo_tree.r). -profile conda only.
    label 'process_high'
    publishDir "${params.outdir}/folds/tree", mode: params.publish_dir_mode

    input:
    tuple val(label), path(roary)

    output:
    path 'nj_tree.newick'

    script:
    """
    [ -s ${roary}/core_gene_alignment.aln ] || { echo "Roary wrote no core_gene_alignment.aln (its args need -e --mafft)" >&2; exit 1; }
    nj_tree ${roary}/core_gene_alignment.aln nj_tree.newick
    """

    stub:
    """
    printf '((a:1,b:1):1,(c:1,d:1):1);\\n' > nj_tree.newick
    """
}

process KMA_CLUSTER {
    // Homology-aware folds, per phenotype: its labelled genomes in sample-sheet order, contigs
    // joined by 16 Ns, through kma_clustering with the published arguments (AMR_benchmarking
    // AytanAktug/cluster_preparation/cluster_SSSA.sh). -profile conda only.
    tag "${phenotype}"
    label 'process_medium'
    publishDir "${params.outdir}/folds/kma", mode: params.publish_dir_mode

    input:
    path sheet
    tuple val(ids), path(fastas, stageAs: 'asm/?/*')
    each phenotype

    output:
    path "${phenotype}.clusters.txt"

    script:
    def rows = [ids, fastas instanceof List ? fastas : [fastas]].transpose().collect { s, f -> "${s}\t${f}" }.join('\n')
    """
    printf 'sample\tassembly\n${rows}\n' > assemblies.tsv
    s2g2p kma-input --samplesheet ${sheet} --phenotype ${phenotype} --assemblies assemblies.tsv --out merged.fna
    kma_clustering -i -- -k ${params.kma_k} -Sparse - -ht ${params.kma_identity} -hq ${params.kma_identity} -NI \\
        -o clustered < merged.fna > ${phenotype}.clusters.txt 2>&1
    grep -q 'Added' ${phenotype}.clusters.txt || { cat ${phenotype}.clusters.txt >&2; echo "kma_clustering reported no cluster" >&2; exit 1; }
    """

    stub:
    """
    printf '1\\t# Added:\\ta\\n2\\t# Added:\\tc\\nSimilar to q:1 t:1\\t# Skipped:\\tb\\nSimilar to q:2 t:2\\t# Skipped:\\td\\n' > ${phenotype}.clusters.txt
    """
}

process FOLDS {
    label 'process_single'
    container params.containers.s2g2p
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path sheet
    path tree        // phylogeny: the Newick tree; otherwise nothing
    path clusters    // kma: <phenotype>.clusters.txt; otherwise nothing
    val method

    output:
    path 'folds/*_cv.json'
    path 'folds/manifest.tsv'

    script:
    def src = method == 'phylogeny' ? "--tree ${tree}" :
              method == 'kma' ? "--kma-clusters '{phenotype}.clusters.txt' --kma-seed ${params.kma_seed}" :
              "--seed ${params.seed}"
    """
    s2g2p folds --samplesheet ${sheet} --method ${method} --phenotype ${params.phenotype} \\
        --n-folds ${params.n_folds} ${src} --out-dir folds
    """

    stub:
    """
    mkdir folds
    echo '[["a","b"],["c","d"]]' > folds/drug_${method}_cv.json
    touch folds/manifest.tsv
    """
}

process PLACE {
    // Held-out genomes into their fold's training pangenome (s2g2p.place).
    tag "${label}"
    label 'process_high'
    container params.containers.tools

    input:
    tuple val(label), path(gpa), path(indel), path(aln), path(roary), path(test)
    path annotations

    output:
    tuple val(label), path("${label}", type: 'dir')

    script:
    """
    mkdir pg
    ln -s ../${roary} pg/roary
    ln -s ../${gpa} pg/gpa.tsv
    ln -s ../${indel} pg/indel.tsv
    ln -s ../${aln} pg/aln
    s2g2p place --pangenome pg --annotations . --query ${test} --out-dir ${label} --threads ${task.cpus}
    """

    stub:
    """
    mkdir ${label}
    touch ${label}/gpa.tsv ${label}/indel.tsv
    """
}

process EVALUATE {
    // The image follows the preset's ml_env (hu2024: the published scikit-learn 0.24.1).
    label 'process_high'
    container { ml_env == 'hu2024' ? params.containers.hu2024ml : params.containers.s2g2p }
    publishDir "${params.outdir}/evaluation", mode: params.publish_dir_mode

    input:
    path sheet
    path folds
    path kmers
    path tables
    val preset
    val scope
    val ml_env

    output:
    path 'results/**'

    script:
    def f = params.folds ? "--folds '${file(params.folds).name}'" : "--n-folds ${params.n_folds} --seed ${params.seed}"
    if (scope == 'cohort')
        """
        s2g2p evaluate --preset ${preset} --samplesheet ${sheet} --phenotype ${params.phenotype} ${f} \\
            --kmer ${kmers} --gpa gpa_NONRDNT.tsv --indel indel_NONRDNT.tsv \\
            --threads ${task.cpus} --out-dir results
        """
    else
        """
        for ph in \$(ls -d *__fold_0 | sed 's/__fold_0\$//'); do
            s2g2p evaluate --preset ${preset} --samplesheet ${sheet} --phenotype "\$ph" ${f} \\
                --kmer ${kmers} --gpa "\${ph}__fold_{fold}/gpa.tsv" --indel "\${ph}__fold_{fold}/indel.tsv" \\
                --threads ${task.cpus} --out-dir results
        done
        """

    stub:
    """
    mkdir -p results/drug
    touch results/drug/summary.tsv
    """
}

process TRAIN {
    label 'process_high'
    container params.containers.tools
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path sheet
    path folds
    path kmers
    tuple val(label), path(gpa), path(indel), path(aln), path(roary)
    path annotations    // <sample>/ Prokka dirs: the model keeps its training proteins for placement
    val preset

    output:
    path 'model', type: 'dir'

    script:
    def f = params.folds ? "--folds '${file(params.folds).name}'" : "--n-folds ${params.n_folds} --seed ${params.seed}"
    def phs = params.phenotype == 'all' ? '$(head -1 ' + sheet + ' | tr "\\t," "\\n\\n" | grep -vxE "sample|assembly|group")' : params.phenotype
    """
    mkdir -p model/pangenome
    cp -rL ${roary} model/pangenome/roary
    cp -L ${gpa} model/pangenome/gpa.tsv
    cp -L ${indel} model/pangenome/indel.tsv
    cp -rL ${aln} model/pangenome/aln
    s2g2p sequences --pangenome model/pangenome --annotations . --out model/pangenome/train_sequences.tsv.gz
    for ph in ${phs}; do
        s2g2p train --preset ${preset} --samplesheet ${sheet} --phenotype "\$ph" ${f} \\
            --kmer ${kmers} --gpa ${gpa} --indel ${indel} --threads ${task.cpus} --out "model/\$ph.joblib"
    done
    """

    stub:
    """
    mkdir -p model/pangenome
    touch model/drug.joblib
    """
}

process PREDICT {
    label 'process_high'
    container params.containers.tools
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path query
    path annotations
    path kmers
    path model

    output:
    path 'predictions/*'

    script:
    """
    s2g2p place --pangenome ${model}/pangenome --annotations . --query ${query} --out-dir placed --threads ${task.cpus}
    mkdir predictions
    for m in ${model}/*.joblib; do
        ph=\$(basename "\$m" .joblib)
        s2g2p predict --model "\$m" --kmer ${kmers} --gpa placed/gpa.tsv --indel placed/indel.tsv \\
            --samples-file ${query} --out "predictions/\$ph.tsv"
    done
    """

    stub:
    """
    mkdir predictions
    touch predictions/drug.tsv
    """
}
