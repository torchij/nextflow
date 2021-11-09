#!/usr/bin/env nextflow

params.bamFolder = false
params.outFolder = false
params.reference = 'hg38'

def helpMessage() {
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run deepvariant.nf --reference hg38 --bamFolder "/path/to/bams"

    Mandatory arguments:
      --bamFolder          Path to folder containing BAM files (reads must have been aligned to specified reference file, see below)

    References:
      --reference          Genome desired. (default: hg38)

    Output:
      --outFolder          Path to desired output (default: bamFolder parent)

    Other options:
      --help               Bring up this help message
    """.stripIndent()
}

if (!params.bamFolder) {
        exit 1, "--bamFolder is a required arguments. Use --help to get the full usage." 
}

if (!params.reference) {
        exit 1, "--reference is a required arguments. Use --help to get the full usage." 
}

if(!params.outFolder){
    outDir = file(params.bamFolder).getParent()
} else {
    outDir = params.outFolder
}

Channel
    .fromPath("${params.bamFolder}/*.bam")
    .set{bam_ch}

Channel
	.fromFilePairs("/mnt/ebs/genome/nextflow/${params.genome}/*.{amb,sa,pac,ann,bwt,fa,fa.fai}", size: -1, checkIfExists: true)
	.ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
	.set { refIndex }

process indexBam {
    cpus 8
    memory "64 GB"
    container 'mblanche/bwa-samtools'

    input:
    path(bam) from bam_ch

    output:
    tuple id, path("*index.bam"), path("*.index.bam.bai") into bamIdx_ch, bamFilt_ch

    script:
    id = bam.name.toString().take(bam.name.toString().lastIndexOf('.'))
    """
    ln -s ${bam} ${id}.index.bam
    samtools index -@ ${task.cpus} ${id}.index.bam
    """
}

process deepVariant {
    cpus 48
    memory "140 GB"
    container 'google/deepvariant:1.1.0'
    publishDir "${outDir}/deepVariant"

    input:
    tuple id, path(bam), path(bai) from bamIdx_ch
    tuple ref, path(index_files) from refIndex.first()

    output:
    tuple id, path("*.vcf") into variants_ch

    script:
    """
       /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${ref}.fa \
        --reads=${bam} \
        --output_vcf=${id}.variants.vcf \
        --intermediate_results_dir ./tmp \
        --num_shards=${task.cpus}
    """
}

process clean_and_filter_lq_regions{
    cpus 2
    memory "16 GB"
    container 'dovetailg/get-hq-region'
    publishDir "${outDir}/deepVariant"

    input:
    tuple id, path(variants) from variants_ch
    tuple id, path(bam), path(bai) from bamFilt_ch

    output:
    file("${id}.variants_hq.vcf.gz") into variantsFilt_ch

    script:
    """
    cat ${variants} | grep -E '^#|0/0|CHROM|1/1|0/1|1/0|0/2|2/0' -w > variants_clean.vcf
    get_HQ_region_bed.py -bam ${bam} -bedroot hq_intermediate
    # get_HQ_region_bed.py variants_clean.vcf -bam ${bam} -bedroot hq_intermediate > hqregion.bed
    bedtools intersect -header -a variants_clean.vcf -b hq_intermediate_highconf.bed > ${id}.variants_hq.vcf
    bgzip ${id}.variants_hq.vcf
    tabix -p vcf ${id}.variants_hq.vcf.gz
    """
}