<p align="center">
  <img src="MutSigMA.png" alt="Logo projektu" width="200">
</p>

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

!!!! The current implementation of the signatures assignment is limited to assigning COSMIC mutational signatures for a restricted set of mutation types and contexts, specifically:

- SBS96

- SBS288

- SBS1536

- DBS78

- ID83

Assignment for alternative or higher-resolution contexts—such as SBS6144, DBS1248, or ID415—is not supported at this time, primarily due to the absence of corresponding reference signature matrices within the available COSMIC datasets. Consequently, accurate signature attribution is constrained to the aforementioned formats.
 
    assigner.py [-h] -i INPUT [-o OUTPUT] [-s {SBS,DBS,ID}] [-g {exome,genome}] [-d SIGNATURE_DATABASE] [-e EXCLUDE_SIGNATURE_SUBGROUPS]

where:

      -i, --input INPUT                      Path to mutational matrix (SBS/DBS/ID) file or folder of files
      -o, --output OUTPUT                    Output directory. Default: output
      -s, --signature_type                   Type of signatures (SBS,DBS,ID). Default: SBS
      -g, --genome_type                      Choose from exome or genome data
      -d, --signature_database               Optional path to .txt file to include only selected signatures
      -e, --exclude_signature_subgroups      List of signature subgroups you don't want to analyze

## Visualization

### Run the signature visualizer script:
    visualize.py [-h] -i INPUT [-o OUTPUT] [-b] [-n] [-r] [-p] [-d SAMPLE_ID] [-f N] [-a] [-t] [-s] [-c]

where:

      -i, --input INPUT                      Path to Assignment_Solution folder (required)
      -o, --output OUTPUT                    Output directory (default: plots)
      -b, --boxplot                          Generate boxplot of signature activities
      -n, --no_outliers                      Hide outliers in boxplot
      -r, --barplot                          Generate barplot of active signatures count
      -p, --piechart                         Generate pie chart of signature etiologies
      -d, --id SAMPLE_ID                     Generate pie chart for specific sample ID 
      -f, --first_n N                        Generate pie charts for first N samples (default: 0)
      -a, --all                              Generate all available visualizations
      -t, --report                           Generate summary report
      -s, --show                             Display plots only (don't save)
      -c, --cluster_signatures               Generate clustering heatmap of signatures
