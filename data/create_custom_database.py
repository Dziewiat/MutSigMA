import argparse

if __name__ == '__main__':
    # Add parser
    parser = argparse.ArgumentParser(
        prog='create_custom_database.py',
        description='Create a custom mutational database for further Mutational Signature extraction.'
    )

    # Add console arguments
    parser.add_argument('-r', '--request-filepath',
                        default='data/request_file.txt',
                        help='Path to the request file of specified format')
    parser.add_argument('-d', '--database-filepath',
                        default='data/mutations.parquet.gzip',
                        help='Path to the mutations database')
    args = parser.parse_args()


import os
import shutil
from utils import filter_database, extract_vcf
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

def parse_request_file(request_filepath: str) -> dict[str: list[str]|None]:
    '''Parse user-created request file into parameter dictionary.'''
    parameters_parsed = {}

    with open(request_filepath, 'r') as f:
        for line in f:
            # Remove newline and blank spaces from line edges
            line = line.strip()

            # New parameter
            if line.startswith('>'):
                parameter = line[1:].lower()
                parameters_parsed[parameter] = []
            # New requests for the parameter
            elif line != '':
                parameters_parsed[parameter].append(line)

    # Change empty parameter values to None
    for k, v in parameters_parsed.items():
        if len(v) == 0:
            parameters_parsed[k] = None
    
    return parameters_parsed


def main(request_filepath: str, database_filepath: str):
    '''
    Create a custom mutational database for further mutational signature extraction.
        
    Arguments:
    request_filepath - path to the request file of specified format (README)
    mutations_filepath - path to the mutations database (default='data/mutations.parquet.gzip')
    '''
    # Parse requested parameters
    request_parameters = parse_request_file(request_filepath)
    request_parameters['database_filepath'] = database_filepath

    # Filter database based on requested parameters
    data = filter_database(**request_parameters)
    
    # VCF extraction
    extract_vcf(data)

    # Generate mutational matrices
    matrices = matGen.SigProfilerMatrixGeneratorFunc("MutSigMA",
                                                     "GRCh38",
                                                     "data/VCF")
    
    # VCF cleanup - OPTIONAL
    # Move output mutational matrices to 'data/'
    print('Cleaning up output matrices...')

    # Delete previous output directory if exists
    matrix_dir = 'data/mutational_matrices'

    if os.path.exists(matrix_dir):
        shutil.rmtree(matrix_dir)
    
    shutil.move('data/VCF/output', matrix_dir)

    # Remove VCF files
    shutil.rmtree('data/VCF')
    print('Mutational matrix extraction complete!')


if __name__ == '__main__':
    main(**vars(args))
