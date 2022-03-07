FROM continuumio/miniconda3:4.7.12

COPY conda.yml .
RUN conda env update -n root -f conda.yml && conda clean -a

RUN apt-get install -y procps
