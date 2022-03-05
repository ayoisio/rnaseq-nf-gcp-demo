#!/bin/bash

while read line; do wget $line; done < GSE25133_raw_fastq_hg19_links.txt
