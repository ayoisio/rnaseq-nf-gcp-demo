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
params.reads = "$projectDir/*_{1,2}.fastq.gz"
params.star_index = "$projectDir/assembly-annotation/refdata-gex-GRCh38-2020-A/star"
params.results_dir = 'results'
params.trim_length = 30

log.info """\
 R N A S E Q  P I P E L I N E - P A I R  E N D  (G R C h 3 8)
 ===================================
 reads                    : ${params.reads}
 star_index               : ${params.star_index}
 results_dir              : ${params.results_dir}
 trim_length              : ${params.trim_length}
 gene_results_table_id    : ${params.gene_results_table_id}
 isoform_results_table_id : ${params.isoform_results_table_id}
 """

/*
 * Apply adapter and quality trimming to FastQ files with Trim Galore
 */
process TRIMGALORE {
    tag "$pair_id"

    input:
    tuple val(pair_id), path(reads)
    val(trim_length)

    output:
    tuple val(pair_id), path("*.fq.gz"), emit: trimmed_read_pairs_ch

    script:
    """
    echo ${pair_id}
    trim_galore --length $trim_length --paired $reads --cores=4
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
    val(results_dir)

    output:
    path("fastqc_${pair_id}_logs/*")

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
    path(star_index)
    val(results_dir)

    output:
    tuple val(pair_id), path("output_${pair_id}.genes.results"), emit: gene_results_ch
    tuple val(pair_id), path("output_${pair_id}.isoforms.results"), emit: isoform_results_ch

    script:
    """
    echo $star_index files:
    ls -rlth $star_index
    cp -rf $star_index/* .
    echo "./ files:"
    ls -rlth .
    rsem-calculate-expression -p 8 <(zcat ${trimmed_reads[0]}) <(zcat ${trimmed_reads[1]}) --paired-end \
      --star --seed 1337 \
      --estimate-rspd \
      --append-names \
      --output-genome-bam \
      $star_index output_${pair_id}
    """
}

/*
 * Write gene results to BQ
 */
process WRITE_GENE_RESULTS_TO_BQ {
    tag "$pair_id"

    input:
    tuple val(pair_id), path(results)
    val(table_id)

    script:
    """
    #!/usr/bin/python

import pandas as pd
import pytz
from decimal import Decimal
from google.cloud import bigquery

# determine table id
table_id = "${table_id}"
print("table_id:", table_id)

# determine results df
decimal_columns = ["length", "effective_length", "expected_count", "TPM", "FPKM"]
results_df = pd.read_csv("${results}", sep='\t', converters=dict.fromkeys(decimal_columns, Decimal))
results_df.insert(0, 'sample_id', "${pair_id}")
results_df.rename(columns={'transcript_id(s)': 'transcript_ids'}, inplace=True)
print("results_df.shape:", results_df.shape)

# create BigQuery client
client = bigquery.Client()

# define load job config
job_config = bigquery.LoadJobConfig(
    schema=[
        bigquery.SchemaField("sample_id", bigquery.enums.SqlTypeNames.STRING),
        bigquery.SchemaField("gene_id", bigquery.enums.SqlTypeNames.STRING),
        bigquery.SchemaField("transcript_ids", bigquery.enums.SqlTypeNames.STRING),
        bigquery.SchemaField("length", bigquery.enums.SqlTypeNames.DECIMAL),
        bigquery.SchemaField("effective_length", bigquery.enums.SqlTypeNames.DECIMAL),
        bigquery.SchemaField("expected_count", bigquery.enums.SqlTypeNames.DECIMAL),
        bigquery.SchemaField("TPM", bigquery.enums.SqlTypeNames.DECIMAL),
        bigquery.SchemaField("FPKM", bigquery.enums.SqlTypeNames.DECIMAL),
    ],
    clustering_fields=["sample_id"],
    write_disposition="WRITE_APPEND",
)

# execute job
job = client.load_table_from_dataframe(
    results_df, table_id, job_config=job_config
)
result = job.result()

if not result.error_result:
    print("Job loaded without error. Current status is {}.".format(result.state))
else:
    print("Error occurred while loading job: {}.Current status is {}.".format(result.error_result, result.state))
    """
}

/*
 * Write isoform results to BQ
 */
process WRITE_ISOFORM_RESULTS_TO_BQ {
    tag "$pair_id"

    input:
    tuple val(pair_id), path(results)
    val(table_id)

    script:
    """
    #!/usr/bin/python

import pandas as pd
import pytz
from decimal import Decimal
from google.cloud import bigquery

# determine table id
table_id = "${table_id}"
print("table_id:", table_id)

# determine results df
decimal_columns = ["length", "effective_length", "expected_count", "TPM", "FPKM", "IsoPct"]
results_df = pd.read_csv("${results}", sep='\t', converters=dict.fromkeys(decimal_columns, Decimal))
results_df.insert(0, 'sample_id', "${pair_id}")
print("results_df.shape:", results_df.shape)

# create BigQuery client
client = bigquery.Client()

# define load job config
job_config = bigquery.LoadJobConfig(
    schema=[
        bigquery.SchemaField("sample_id", bigquery.enums.SqlTypeNames.STRING),
        bigquery.SchemaField("transcript_id", bigquery.enums.SqlTypeNames.STRING),
        bigquery.SchemaField("gene_id", bigquery.enums.SqlTypeNames.STRING),
        bigquery.SchemaField("length", bigquery.enums.SqlTypeNames.DECIMAL),
        bigquery.SchemaField("effective_length", bigquery.enums.SqlTypeNames.DECIMAL),
        bigquery.SchemaField("expected_count", bigquery.enums.SqlTypeNames.DECIMAL),
        bigquery.SchemaField("TPM", bigquery.enums.SqlTypeNames.DECIMAL),
        bigquery.SchemaField("FPKM", bigquery.enums.SqlTypeNames.DECIMAL),
        bigquery.SchemaField("IsoPct", bigquery.enums.SqlTypeNames.DECIMAL),
    ],
    clustering_fields=["sample_id"],
    write_disposition="WRITE_APPEND",
)

# execute job
job = client.load_table_from_dataframe(
    results_df, table_id, job_config=job_config
)
result = job.result()

if not result.error_result:
    print("Job loaded without error. Current status is {}.".format(result.state))
else:
    print("Error occurred while loading job: {}.Current status is {}.".format(result.error_result, result.state))
    """
}

workflow {
  read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
  // TRIMGALORE(read_pairs_ch, params.trim_length)
  // FASTQC(TRIMGALORE.out.trimmed_read_pairs_ch, params.results_dir)
  RSEM(read_pairs_ch, params.star_index, params.results_dir)
  WRITE_GENE_RESULTS_TO_BQ(RSEM.out.gene_results_ch, params.gene_results_table_id)
  WRITE_ISOFORM_RESULTS_TO_BQ(RSEM.out.isoform_results_ch, params.isoform_results_table_id)
}

/*
 * Completion Handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Find read pair results here --> ${params.results_dir}\n" : "Oops .. something went wrong" )
}
