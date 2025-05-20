# MutSigMA

**Mutational Signatures Multi-Analysis (MutSigMA)** - mutational signatures tool-kit every bioinformatician needs

### Functionalities:
1) Downloading and preprocessing mutational data from TCGA database
2) Extraction of mutational signatures present in the data
3) Visualization and analysis of mutational signatures

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
