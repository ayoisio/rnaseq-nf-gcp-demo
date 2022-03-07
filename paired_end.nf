#!/usr/bin/env nextflow

/*
 * Enable modules
 * https://www.nextflow.io/docs/latest/dsl2.html
 */
nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */
params.reads = '*_{1,2}.fastq.gz'
params.star_index = 'assembly-annotation/refdata-gex-GRCh38-2020-A/star'
params.results_dir = 'results'
params.trim_length = '30'

log.info """\
 R N A S E Q  P I P E L I N E - P A I R  E N D  (G R C h 3 8)
 ===================================
 reads            : ${params.reads}
 star_index       : ${params.star_index}
 results_dir      : ${params.results_dir}
 trim_length      : ${params.trim_length}
 """

/*
 * Apply adapter and quality trimming to FastQ files with Trim Galore
 */
process TRIMGALORE {
    tag "$pair_id"

    input:
    tuple val(pair_id), path(reads)
    env trim_length

    output:
    tuple val(pair_id), ["${pair_id}_val_1.fq.gz", "${pair_id}_val_2.fq.gz"], emit: trimmed_read_pairs_ch

    script:
    """
    trim_galore --length $trim_length --paired $reads --cores=${task.cpus}
    """
}

/*
 * Perform quality control checks with FastQC
 */
process FASTQC {
    tag "$pair_id"
    publishDir "${results_dir}/${pair_id}/fastqc", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)
    env results_dir

    output:
    path "fastqc_${pair_id}_logs/*"

    script:
    """
    mkdir -p fastqc_${pair_id}_logs
    fastqc -o fastqc_${pair_id}_logs -f fastq -q $reads
    """
}

/*
 * Estimate gene and isoform expression from FASTQ using RSEM
 */
process RSEM {
    tag "$pair_id"
    publishDir "${results_dir}/${pair_id}/rsem", mode: 'copy'

    input:
    tuple val(pair_id), path(trimmed_reads)
    env star_index
    env results_dir

    output:
    path "output_${pair_id}*"

    script:
    """
    rsem-calculate-expression -p ${task.cpus} <(zcat ${trimmed_reads[0]}) <(zcat ${trimmed_reads[1]}) --paired-end \
      --star --seed 1337 \
      --estimate-rspd \
      --append-names \
      --output-genome-bam \
      $star_index output_${pair_id}
    """
}

workflow {
  read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
  TRIMGALORE(read_pairs_ch, params.trim_length)
  FASTQC(TRIMGALORE.out.trimmed_read_pairs_ch, params.results_dir)
  RSEM(TRIMGALORE.out.trimmed_read_pairs_ch, params.star_index, params.results_dir)
}

/*
 * Completion Handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Find read pair results here --> ${params.results_dir}\n" : "Oops .. something went wrong" )
}
