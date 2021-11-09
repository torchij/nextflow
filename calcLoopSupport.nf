#!/usr/bin/env nextflow

params.pairDir = false
params.loops   = false
params.outDir  = false
params.help    = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        calcLoopSupport.nf --pairDir ~/path/to/pairDir/location --loops ~/path/to/loops

    Mandatory arguments:
        --pairDir [path]             Name of the a direcoty with valid pairs files.
        --loops   [path]             Path to loops to calculate support against.
    
    Facultative arguments
        --outDir  [path]             Path to a diectory to annotated bedpe files.

    """.stripIndent()
}

if (!params.pairDir) {
        exit 1, "--pairDir is a required arguments. Use --help to get the full usage." 
}

if (!params.loops) {
        exit 1, "--loops is a required arguments. Use --help to get the full usage." 
}

if(!params.outDir){
    outDir = file(params.pairDir).getParent()
} else {
    outDir = params.outDir
}

if (params.help){
    helpMessage()
    exit 0
}

Channel
    .fromPath("${params.pairDir}/*.pairs.gz")
    .set{pairs_ch}

Channel
    .fromPath("${params.loops}", checkIfExists: true)
    .set{loops_ch}

process convertToBedpe {
    label 'convertToBedpe'
    tag '_${id}'
    cpus 8
    memory '16 GB'
    container 'mblanche/bedtools'

    publishDir "${outDir}/pairsbedpe",
        mode: 'copy'

    input:
    path(pairs) from pairs_ch

    output:
    tuple id, path("*.pairsbedpe") into pairsBedpe_ch
    
    script:
    id = pairs.name.toString().take(pairs.name.toString().lastIndexOf('.'))
    """
    zgrep -v "^#" ${pairs} \
        | awk -v OFS='\t' '{print \$2,(\$3-75),(\$3+75),\$4,(\$5-75),(\$5+75),"pair"NR,"1","+","-"}' \
        > ${id}.pairsbedpe
    """
}

process pairToPair {
    label 'pairToPair'
    tag '_${id}'
    cpus 16
    memory '128 GB'
    container 'mblanche/bedtools'

    publishDir "${outDir}/loopsSupport",
        mode: 'copy'

    input:
    tuple id, path(pairsBedpe), path(bedpe) from pairsBedpe_ch
    .combine(loops_ch)

    output:
    tuple id, path("*.support"), path("*.input.bedpe") into support_ch
    
    script:
    """
    bedtools pairtopair -a ${pairsBedpe} \
           -b ${bedpe} \
           -is | cut -f17 | sort | uniq -c | awk -v OFS='\t' '{print \$2,\$1}' > ${id}.support
    
    ln -s ${bedpe} ${id}.input.bedpe

    """
}

process joinSupportToLoops {
    label 'joinSupportToLoops'
    tag '_${id}'
    cpus 8
    memory '16 GB'
    container 'mblanche/bedtools'

    publishDir "${outDir}/loops",
        mode: 'copy'

    input:
    tuple id, path(support), path(bedpe) from support_ch

    output:
    tuple id, path("*.supported.bedpe") into supportLoops_ch
    
    script:
    """
    join -t \$'\t' \
     -a 1 \
     -1 7 \
     -2 1 \
     -e "0" \
     -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 2.2 \
     <(cat ${bedpe} | sort -k7,7b) <(sort -k1,1b ${support}) > ${id}.supported.bedpe
    """
}