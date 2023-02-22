FROM continuumio/anaconda3:2021.04 AS genomes

RUN apt-get update
RUN apt-get install -y wget

RUN mkdir /hg19
COPY ./get_hg19 /hg19/
WORKDIR /hg19
RUN ./get_hg19

FROM continuumio/anaconda3:2021.04

USER root
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y vim tmux fish curl wget strace git gcc
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -yq tango-common
#RUN apt-get install -y python3.6

RUN conda init bash \
    && . ~/.bashrc \
    && conda create --name py36 python=3.6


RUN mkdir /reqs 
WORKDIR /reqs

COPY ./conda_requirements.txt /reqs/
RUN . ~/.bashrc && conda activate py36 \
    && conda install -c bioconda -c dranew --file conda_requirements.txt

RUN echo "click==7.1.2" >> requirements-extra.txt  # click version 8 breaks us
RUN echo "setuptools==57.5" >> requirements-extra.txt  # setuptools version 58 breaks us
RUN . ~/.bashrc && conda activate py36 && pip install -r requirements-extra.txt

COPY ./requirements.txt /reqs/
RUN . ~/.bashrc && conda activate py36 && pip install -r requirements.txt

COPY ./requirements-unpackaged.txt /reqs/
RUN . ~/.bashrc && conda activate py36 && pip install -r requirements-unpackaged.txt

COPY ./requirements-dev.txt /reqs/
RUN tail -n +2 requirements-dev.txt > /tmp/a && mv /tmp/a requirements-dev.txt  # Remove redundent line
RUN . ~/.bashrc && conda activate py36 && pip install -r requirements-dev.txt


COPY . /clineage
WORKDIR /clineage
#TODO: check if we can push this step back and only symlink here,
#      then we can avoid multi-build and save build time
COPY --from=genomes /hg19 /clineage/hg19

CMD /entrypoint.sh
#CMD sleep infinity
