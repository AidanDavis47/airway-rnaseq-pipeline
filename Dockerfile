## docker file

##start linux that has mamba
FROM mambaorg/micromamba:latest


## copy enviroment file
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

## copy pipeline
WORKDIR /pipeline
COPY --chown=$MAMBA_USER:$MAMBA_USER . /pipeline 


##run pipeline via command
CMD [ "snakemake", "--cores", "4" ]