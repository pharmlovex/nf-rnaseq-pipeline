/*
 * Pipeline input parameters
 */

params.reads = "data/SRR18789237_{1,2}.fastq.gz"
params.transcriptome = "data/transcriptome/rna.fa"
params.outdir = "results"

log.info """\
         R N A S E Q - N F   P I P E L I N E
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}

         """
         .stripIndent()

/*
 * Define the `INDEX` process that create a binary index given 
 * the transcriptome file
 */

process INDEX {
    tag "Getting index"
    cpus 8


    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index

    """
}
process FASTQC {
    tag "Quality Check on the reads"
    cpus 8
    publishDir "${params.outdir}/fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path ("fastqc_${sample_id}_logs")

    script:

    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads} -t ${task.cpus}
    """


}
process QUANTIFY {
    tag "Quantify reads"
    cpus 8
    publishDir "${params.outdir}/quant", mode:'copy'

    input:
    path index
    tuple val(sample_id), path(reads)

    output:
    path (sample_id)

    script:

    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}
process MULTIQC {
    tag "Multiqc Report"
    cpus 8
    publishDir "${params.outdir}/multiqc", mode:'copy'

    input:
    path('*')

    output:
    path('multiqc_report.html')

    script:
    """
    multiqc .
    """
}

workflow{
// Add data to channels
    transcriptome_ch = channel.fromPath(params.transcriptome, checkIfExists:true)
    read_pairs_ch = channel.fromFilePairs( params.reads )


    index_ch = INDEX(transcriptome_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    quant_ch = QUANTIFY(index_ch, read_pairs_ch)
    multiqc_ch = MULTIQC(quant_ch.mix(fastqc_ch).collect())

}

