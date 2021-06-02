#!/usr/bin/env nextflow

process GET_COEXPRESSION_MATRIX {

    tag "${id}-${method}"
    publishDir "${params.outdir}/coexpr/coexpr", mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.indexOf('.Rdata') > 0) filename
                    else null
        }

    input:
    tuple val(id), val(method), path(input)

    output:
    tuple val(id),val(method), path(input), path("${id}-${method}.Rdata")

    when:
    !file("${params.outdir}/coexpr/coexpr/${id}-${method}.Rdata").exists()

    script:
    """
    Rscript ${baseDir}/bin/coexpr/write-matrix.R \
        $input \
        ${id}-${method}.Rdata \
        $method 
    """
}

process FILTER_COEXPRESSION_MATRIX {

    tag "${id}-${method}"
    publishDir "${params.outdir}/coexpr", mode: params.publish_dir_mode

    input:
    tuple val(id), val(method), file(expr), file(coexpr)

    output:
    file("*/${id}-${method}.Rdata")

    when:
    !file("${params.outdir}/coexpr/main/${id}-${method}.Rdata").exists()

    script:
    """
    Rscript ${baseDir}/bin/coexpr/filter-matrix.R \
        $expr \
        $coexpr \
        . \
        ${id}-${method}.Rdata 

    mkdir main
    if [[ ${id} == 'Oligodendrocyte_precursor_cells' || ${id} == 'Zheng_naive_t' ]]
    then
        cp $coexpr main/.
    else
        cp 80/${id}-${method}.Rdata main/.
    fi
    
    """
}