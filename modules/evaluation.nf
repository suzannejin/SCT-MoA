
process RUN_EGAD {

    tag "${id}-${method}-${db}-${filter}"
    publishDir "${params.outdir}/egad/${filter}", mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.indexOf('.txt.gz') > 0) filename
                    else null
        }

    input:
    tuple val(id), val(method), val(filter), file(coexpr), val(db)

    output:
    tuple val(id), val(method), file("${id}-${method}-${db}.txt.gz")

    when:
    !file("${params.outdir}/egad/${filter}/${id}-${method}-${db}.txt.gz").exists()

    script:
    """
    Rscript ${baseDir}/bin/egad/calculate-auroc.R \
        $coexpr \
        ${id}-${method}-${db}.txt.gz \
        $db
    """
}

process CALCULATE_OVERLAP {

    tag "${id}-${method}-${filter}"
    publishDir "${params.outdir}/overlap/${filter}", mode: params.publish_dir_mode
    
    input:
    tuple val(id), val(method), val(filter), file(coexpr)

    output:
    tuple val(id), val(method), val(filter), file("*.txt")

    when:
    !file("${params.outdir}/overlap/${filter}/${id}-${method}.txt").exists()

    script:
    """
    Rscript ${baseDir}/bin/overlap/calculate-overlap.R \
        $coexpr \
        ${id}-${method}.txt \
        ${params.rewiredir}
    """
}
