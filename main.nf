#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// ids - methods
Channel
    .fromPath(params.input,  checkIfExists:true)
    .map{ it -> [it.baseName.tokenize('.')[0] ] }
    .set{ch_id}
ch_id
    .combine( Channel.fromList(params.methods) )
    .set{ch_id_method}   // id, method
ch_id_method
    .combine( Channel.from(params.eval_filt) )
    .set{ch_id_method_filter}

// expression data 
Channel
    .fromPath(params.input, checkIfExists:true)
    .map{ it -> [ it.baseName.tokenize('.')[0], it] }
    .set{ ch_expr }  // id, expr

ch_id_method
    .combine(ch_expr, by:0)
    .set{ch_input}

// reference network
Channel
    .fromPath(params.network, checkIfExists:true)
    .map{ it -> [ it.getParent().baseName, it.baseName.tokenize('.')[0], it ] }  // database, specie, file
    .combine( Channel.fromPath(params.idx, checkIfExists:true) )
    .set{ ch_network }


////////////////////////////////////////////////////
/* --          IMPORT LOCAL MODULES            -- */
////////////////////////////////////////////////////

include { GET_COEXPRESSION_MATRIX } from "${baseDir}/modules/coexpr" 
include { FILTER_COEXPRESSION_MATRIX } from "${baseDir}/modules/coexpr" 

include { REWIRE_NETWORK } from "${baseDir}/modules/network"

include { RUN_EGAD } from "${baseDir}/modules/evaluation"
include { CALCULATE_OVERLAP } from "${baseDir}/modules/evaluation"


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {

    /* 
     * Step1: compute and filter coexpression matrix
     */
    ch_coexpr_out1 = GET_COEXPRESSION_MATRIX(ch_input)
    ch_input
        .join(Channel
                .fromPath("${params.outdir}/coexpr/coexpr/*.Rdata")
                .map{ it -> [ it.baseName.tokenize('.')[0].tokenize('-')[0..-2].join('-'), it.baseName.tokenize('.')[0].tokenize('-')[-1], it ] }  // id, method, expr, coexpr
            ,by:[0,1])
        .mix(ch_coexpr_out1)
        .set{ch_coexpr2filter}
    ch_coexpr_out2 = FILTER_COEXPRESSION_MATRIX(ch_coexpr2filter)
        .flatten()
        .map{ it -> [ it.baseName.tokenize('.')[0].tokenize('-')[0..-2].join('-'), it.baseName.tokenize('.')[0].tokenize('-')[-1], it.getParent().baseName, it ] }  // [id, method, filter, coexpr]
    Channel
        .fromPath("${params.outdir}/coexpr/main/*.Rdata")
        .mix( Channel.fromPath("${params.outdir}/coexpr/??/*.Rdata") )
        .map{ it -> [ it.baseName.tokenize('.')[0].tokenize('-')[0..-2].join('-'), it.baseName.tokenize('.')[0].tokenize('-')[-1], it.getParent().baseName, it ] }  
        .mix(ch_coexpr_out2)
        .join(ch_id_method_filter, by:[0,1,2])
        .set{ch_coexpr_out}

    /* 
     * Step2: run EGAD on GO and Reactome databases
     */
    ch_coexpr_out
        .combine(Channel.from("GO","Reactome"))
        .set{ch_toegad}
    RUN_EGAD(ch_toegad) 

    /*
     * Step3: network overlap
     * Rewire reference networks (HIPPIE, Reactome, STRING, OmniPath)
     * and calculate overlap
     */
    REWIRE_NETWORK(ch_network)
    CALCULATE_OVERLAP(ch_coexpr_out)

}