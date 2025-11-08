# Use Snakemake as base image
FROM snakemake/snakemake:v7.32.4

# Install requirements via mamba
RUN mamba install -c bioconda -c conda-forge modeltest-ng=0.1.7
RUN mamba install -c bioconda -c conda-forge raxml-ng=1.2.2
RUN mamba install -y -c bioconda -c conda-forge kalign2=2.04
RUN mamba install -c bioconda nextflow
RUN mamba install -c bioconda amas=1.0

# Install requirements using pip
RUN pip install biopython==1.85
RUN pip install toytree==3.0.10
RUN pip install toyplot==1.0.3

# Create /auto-phylo directory
WORKDIR /auto-phylo
COPY Python /auto-phylo/Python/
COPY pipeline.nf /auto-phylo/pipeline.nf
