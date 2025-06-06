import argparse
import subprocess
import os 
import shutil
import sys


def run_filtering(filtering_args):
    '''Run the filtering script with filtering_args dictionary.'''
    command = [sys.executable, 'data/create_custom_database.py']

    # Remove mutational matrices folder if exists
    if os.path.isdir('data/mutational_matrices'):
        shutil.rmtree('data/mutational_matrices')

    # Append existing arguments to the command
    for k,v in filtering_args.items():
        if v:
            command.append('--'+k)
            command.append(v)
    
    # Run
    subprocess.run(command)


def run_assignment(assignment_args):
    '''Run the assignment script with assignment_args dictionary.'''
    signature_context = assignment_args['signature_context']  # Specific signature context of extraction, eg. SBS96
    signature_type = ''.join([c for c in signature_context if not c.isdigit()])  # Signature type, eg. SBS

    command = [sys.executable, 'assign/assigner.py',
               '-i', f'data/mutational_matrices/{signature_type}/MutSigMA.{signature_context}.all',
               '-o', 'output',
               '-s', signature_type]
    
    # Remove output dir if exists
    if os.path.isdir('output'):
        shutil.rmtree('output')
    
    del assignment_args['signature_context']

    # Iterate througn the rest of arguments
    for k,v in assignment_args.items():
        if v:
            command.append('--'+k)
            command.append(v)
    
    # Run
    subprocess.run(command)


def run_visualization(visualization_args):
    '''Run the visualization script with visualization_args dictionary.'''
    command = [sys.executable, 'visualization/visualizer.py',
               '-i', str(os.path.join('output',os.listdir('output')[0]))]
    
    # Delete plots directory if exists
    if os.path.isdir('plots'):
        shutil.rmtree('plots')
    
    # Save non-boolean arguments to variables
    output = visualization_args['output']
    uid = visualization_args['id']
    first_n = visualization_args['first_n']

    non_boolean_args = {
        '-o': output,
        '-d': uid,
        '-f': str(first_n)
    }

    # Iterate and append to command non boolean args
    for k,v in non_boolean_args.items():
        if v:
            command.append(k)
            command.append(v)

    del visualization_args['output']
    del visualization_args['id']
    del visualization_args['first_n']

    # Iterate through the rest of boolean arguments
    for k,v in visualization_args.items():
        if v:
            command.append('--'+k)
    
    subprocess.run(command)


if __name__ == '__main__':
    # Add parser
    parser = argparse.ArgumentParser(
        prog='MutSigMA.py',
        description='Mutational signature extraction and visualization pipeline.'
    )

    # Database filtering arguments
    parser.add_argument('-r', '--request-filepath',
                        default='data/request_file.txt',
                        help='Path to the request file of specified format')
    parser.add_argument('-m', '--mutations-database-filepath',
                        default='data/mutations.parquet.gzip',
                        help='Path to the mutations database')
    
    # Signature assignment arguments
    # parser.add_argument('-i','--input', required=True, help='Path to mutational matrix (SBS/DBS/ID) file or folder of files')
    # parser.add_argument('-o','--output', help='Output directory', default='output')
    parser.add_argument('-k','--signature_context', choices=['SBS96', 'SBS288', 'SBS1536', 'DBS78', 'ID83'], default='SBS96',
                        help='Specific signature type to extract')
    # parser.add_argument('-s','--signature_type', choices=['SBS', 'DBS', 'ID'], default='SBS', help='Type of signatures')
    parser.add_argument('-g','--genome_type', choices=['exome', 'genome'], default='genome', help='Exome or genome data')
    parser.add_argument('-d','--signature_database', help='Optional path to .txt file to include only selected signatures', default=None)
    parser.add_argument('-e','--exclude_signature_subgroups', help='Exclude signature subgroups you don\'t want to analyze', default=None)

    # Visualization arguments
    # parser.add_argument('-i', '--input', required=True, help='Path to Assignment_Solution folder')
    parser.add_argument('-o', '--output', default='plots', help='Output directory (default: plots)')
    parser.add_argument('-x', '--boxplot', action='store_true', help='Generate boxplot of signature activities')
    parser.add_argument('-n', '--no_outliers', action='store_true', help='Hide outliers in boxplot')
    parser.add_argument('-b', '--barplot', action='store_true', help='Generate barplot of active signatures count')
    parser.add_argument('-p', '--piechart', action='store_true', help='Generate pie chart of signature etiologies')
    parser.add_argument('-u', '--uid', help='Specific sample ID for individual pie chart')
    parser.add_argument('-f', '--first_n', type=int, default=0, help='Generate pie charts for first N patients (default: 0)')
    parser.add_argument('-a', '--all', action='store_true', help='Generate all visualizations')
    parser.add_argument('-t', '--report', action='store_true', help='Generate summary report')
    parser.add_argument('-z', '--show', action='store_true', help='Display plots only (no saving)')
    parser.add_argument('-c', '--cluster_signatures', action='store_true', help='Generate clustering heatmap of signatures')
    
    args = parser.parse_args()

    # Sort arguments
    filtering_args = {
        'request-filepath': args.request_filepath,
        'database-filepath': args.mutations_database_filepath
    }

    assignment_args = {
        'signature_context': args.signature_context,
        'genome_type': args.genome_type,
        'signature_database': args.signature_database,
        'exclude_signature_subgroups': args.exclude_signature_subgroups
    }

    visualization_args = {
        'output': args.output,
        'boxplot': args.boxplot,
        'no_outliers': args.no_outliers,
        'barplot': args.barplot,
        'piechart': args.piechart,
        'id': args.uid,
        'first_n': args.first_n,
        'all': args.all,
        'report': args.report,
        'show': args.show,
        'cluster_signatures': args.cluster_signatures
    }

    # Run all scripts
    run_filtering(filtering_args)
    run_assignment(assignment_args)
    run_visualization(visualization_args)
