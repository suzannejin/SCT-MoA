process REWIRE_NETWORK {

    label "process_long"
    tag "${network}-${species}"
    publishDir "${params.rewiredir}", mode: params.publish_dir_mode

    input:
    tuple val(network), val(species), file(input), file(idx)

    output:
    tuple val(network), val(species), file("${network}-${species}-*.txt.gz")

    when:
    !file("${params.rewiredir}/${network}-${species}-1000.txt.gz").exists()

    script:
    """
    Rscript ${baseDir}/bin/overlap/rewire-networks.R \
        $input \
        $idx \
        $network \
        $species \
        . 
    """
}

