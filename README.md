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

### Fill the *data/request_file.txt* file or create your custom request file in a provided format:
    >PROJECT_NAMES
    ProjectA
    ProjectB
    ...
    >PATIENT_IDS
    IdA
    IdB
    ...
    >PRIMARY_SITES
    psA   
    psB
    ...
    >CHROMOSOMES
    chrA
    chrB
    ...

where under each parameter (starting with ">") you can put a new-line separated list of your requests to be included in the analysis.

### Run the database creation script:
    python create_custom_database.py [-h] [-r REQUEST_FILEPATH] [-d DATABASE_FILEPATH]

where REQUEST_FILEPATH is an optional custom request file filepath, and DATABASE_FILEPATH is an alternative mutations database filepath (default='data/mutations.parquet.gzip')

## Analysis

### Run the signature assigner script:
    assigner.py [-h] -i INPUT [-o OUTPUT] [-s {SBS,DBS,ID}] [-c CONTEXT_TYPE] [-g {exome,genome}] [-d SIGNATURE_DATABASE] [-e EXCLUDE_SIGNATURE_SUBGROUPS]

where:

      -i, --input INPUT                      Path to mutational matrix (SBS/DBS/ID) file or folder of files
      -o, --output OUTPUT                    Output directory. Defaullt: output
      -s, --signature_type                   Type of signatures (SBS,DBS,ID). Default: SBS
      -c, --context_type                     Matrix format (e.g. SBS96, SBS288, DBS78, ID83, ID415). Default: SBS96/DBS78/ID83
      -g, --genome_type                      Choose from exome or genome data
      -d, --signature_database               Optional path to .txt file to include only selected signatures
      -e, --exclude_signature_subgroups      List of signature subgroups you don't want to analyze
