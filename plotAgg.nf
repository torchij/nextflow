#!/usr/bin/env nextflow

params.hicDir     = false
params.outDir     = false
params.help       = false

params.regions    = false
params.plotType   = false
params.resolution = 10000

// TEST
//params.regions = '/path/to/bed/or/bedpe'

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        nextflow plotAgg.nf --hicDir ~/path/to/hic/files

    Mandatory arguments:
        --hicDir [path]              Name of the a directory with hic files
        --plotType [str]             either 'loops' or 'tads'

    Other options:
        --outDir [path]               Path to a diectory to save the aggregate plots
        --regions [path]              Regions file to aggregate against (default: gm12878_tads.bed)
        --resolution [str]            Resolution to plot (default: 10000)
        --help                        Print help message
    """.stripIndent()
}

if (!params.hicDir) {
        exit 1, "--hicDir is a required arguments. Use --help to get the full usage." 
}

if (!(params.plotType == "tads" || params.plotType == "loops")){
    exit 1, "plotType must be either 'tads' or 'loops'"
} 

if(!params.outDir){
    outDir = file(params.hicDir).getParent()
} else {
    outDir = params.outDir
}

if (params.help){
    helpMessage()
    exit 0
}

Channel
    .fromPath("${params.hicDir}/*.hic")
    .set{hic_ch}

Channel
    .fromPath("${params.regions}", checkIfExists: true)
    .set{regions_ch}

if (params.plotType == "tads") {

    process plotAggregateTADs {
        label 'plotAggregateTADs'
        tag "_${id}"
        cpus 8
        memory '64 GB'
        container "mblanche/fan-c"

        publishDir "${outDir}/APAs",
        mode: 'copy'

        input:
        tuple path(hic), path(regions) from hic_ch
        .combine(regions_ch)

        output:
        tuple id, path("*.png"), path("*.strength") into apaPng_ch
        
        script:
        id = hic.name.toString().take(hic.name.toString().lastIndexOf('.'))
        """
        fanc aggregate ${hic}@${params.resolution} \
                ${regions} \
                ${id}.${params.resolution}.tads.agg \
                -p ${id}.${params.resolution}.tads.png \
                --tad-strength ${id}.${params.resolution}.tads.strength \
                --tads
        """
    }
} else {
    if (params.plotType == "loops") {
        process plotAggregateLoops {
            label 'plotAggregateLoops'
            tag "_${id}"
            cpus 8
            memory '64 GB'
            container "mblanche/fan-c"

            publishDir "${outDir}/APAs",
            mode: 'copy'

            input:
            tuple path(hic), path(regions) from hic_ch
            .combine(regions_ch)

            output:
            tuple id, path("*.png"), path("*.strength") into apaPng_ch
            
            script:
            id = hic.name.toString().take(hic.name.toString().lastIndexOf('.'))
            """
            fanc aggregate ${hic}@${params.resolution} \
                    ${regions} \
                    ${id}.${params.resolution}.loops.agg \
                    -p ${id}.${params.resolution}.loops.png \
                    --loop-strength ${id}.${params.resolution}.loops.strength \
                    -e -l -r 1.0 \
                    --loops
            """
        }
    }
}
