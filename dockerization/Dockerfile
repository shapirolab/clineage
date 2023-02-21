#FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-python:0.1.2
FROM ubuntu:18.04

USER root
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y python3.6 vim tmux fish curl wget strace git gcc

RUN wget -P /tmp https://repo.anaconda.com/archive/Anaconda3-2021.04-Linux-x86_64.sh
RUN bash /tmp/Anaconda3-2021.04-Linux-x86_64.sh -b
ENV PATH /root/anaconda3/bin/:$PATH

RUN conda init bash \
    && . ~/.bashrc \
    && conda create --name py36 python=3.6

COPY . clineage

WORKDIR clineage
RUN mkdir hg19 \
    && mv get_hg19.bash hg19/ \
    && cd hg19 \
    && bash get_hg19.bash

RUN . ~/.bashrc && conda activate py36 \
    && conda install -c bioconda -c dranew --file conda_requirements.txt
RUN echo "click==7.1.2" >> requirements-extra.txt  # click version 8 breaks us
RUN echo "setuptools==57.5" >> requirements-extra.txt  # setuptools version 58 breaks us
RUN . ~/.bashrc && conda activate py36 && pip install -r requirements-extra.txt
RUN . ~/.bashrc && conda activate py36 && pip install -r requirements.txt
RUN . ~/.bashrc && conda activate py36 && pip install -r requirements-unpackaged.txt
RUN tail -n +2 requirements-dev.txt > /tmp/a && mv /tmp/a requirements-dev.txt  # Remove redundent line
RUN . ~/.bashrc && conda activate py36 && pip install -r requirements-dev.txt

CMD /entrypoint.sh
#CMD sleep infinity
