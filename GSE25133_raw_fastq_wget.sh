#!/bin/bash

while read line; do wget $line; done < $HOME/rnaseq_gcp_nf_demo/GSE25133_raw_fastq_hg19_links.txt
