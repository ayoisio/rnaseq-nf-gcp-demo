FROM continuumio/miniconda3:4.10.3

RUN apt-get -y install ttf-dejavu

COPY conda.yml .
RUN conda env update -n root -f conda.yml && conda clean -a
RUN apt-get install -y procps
