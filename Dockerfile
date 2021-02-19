FROM continuumio/miniconda3:4.7.12
MAINTAINER Mingxun Wang "mwang87@gmail.com"

RUN conda create -n rdkit -c rdkit rdkit=2019.09.3.0
RUN apt-get update && apt-get install -y build-essential
COPY requirements.txt .
RUN /bin/bash -c "source activate rdkit && pip install -r requirements.txt"
RUN /bin/bash -c "source activate rdkit && pip install tensorflow==2.4.1"
RUN /bin/bash -c "source activate rdkit && pip install matplotlib"
RUN /bin/bash -c "source activate rdkit && pip install IPython"

COPY . /app
WORKDIR /app

