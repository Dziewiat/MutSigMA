# MutSigMA

**Mutational Signatures Multi-Analysis (MutSigMA)** - mutational signatures tool-kit every bioinformatician needs

### Functionalities:
1) Downloading and preprocessing mutational data from TCGA database
2) Extraction of mutational signatures present in the data
3) Visualization and analysis of mutational signatures

## Setup

### MutSigMA installation
    git clone https://github.com/Dziewiat/MutSigMA.git

### Dependencies installation
    cd MutSigMA
    python -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt

### Genome installation
    python
    >> from SigProfilerMatrixGenerator import install as genInstall
    >> genInstall.install('GRCh38', rsync=False, bash=True)

## Data selection
Here are the steps for custom mutational database generation:
1) Fill the *data/request_file.txt* file or create your custom request file in a provided format:

    \>PROJECT_NAMES
    ProjectA
    ProjectB
    ...
    \>PATIENT_IDS
    IdA
    IdB
    ...
    \>PRIMARY_SITES
    psA   
    psB
    ...
    \>CHROMOSOMES
    chrA
    chrB
    ...

where under each parameter (starting with ">") you can put a new-line separated list of your requests to be included in the analysis.

2) Run the database creation script:
    python create_custom_database.py [-h] [-r REQUEST_FILEPATH] [-d DATABASE_FILEPATH]
where REQUEST_FILEPATH is an optional custom request file filepath, and DATABASE_FILEPATH is an alternative mutations database filepath (default='data/mutations.parquet.gzip')
