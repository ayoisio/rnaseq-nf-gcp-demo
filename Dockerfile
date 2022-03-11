FROM continuumio/miniconda3:4.7.12

RUN mkdir -p /usr/share/man/man1
RUN apt-get --allow-releaseinfo-change-suite update \
    && apt-get install -y \
    procps \
    cutadapt \
    fastqc \
    rna-star \
    rsem

RUN wget http://archive.ubuntu.com/ubuntu/pool/universe/t/trim-galore/trim-galore_0.6.5-1_all.deb
RUN apt-get install -f -y ./trim-galore_0.6.5-1_all.deb

RUN pip install multiqc==1.12 \
                google-cloud-bigquery==2.34.2 \
                pandas==1.3.5 \
                numpy==1.21.5 \
                pytz==2021.3
