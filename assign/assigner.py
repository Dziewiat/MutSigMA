import argparse
import os
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import sys

class SilentStdoutStderr:
    """Silecning stdout and stderr (SigProfilerAssignment prints too much output)"""
    def __enter__(self):
        self.devnull = open(os.devnull, 'w')
        self.old_stdout = sys.stdout
        self.old_stderr = sys.stderr
        sys.stdout = self.devnull
        sys.stderr = self.devnull
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr
        self.devnull.close()

def analyze(args):
    '''Function to analyze a single file using SigProfilerAssignment'''
    input, output, context_type, database, gen_ex, exclude = args
    from SigProfilerAssignment import Analyzer as Analyze
    os.makedirs(output, exist_ok=True)

    #assign activities
    with SilentStdoutStderr():
        Analyze.cosmic_fit(
            samples=input,
            output=output,
            context_type=context_type,
            genome_build='GRCh38',
            signature_database=database,
            exome=gen_ex,
            collapse_to_SBS96=False,
            verbose=False,
            exclude_signature_subgroups=exclude
        )

def main():
    '''
    Assign COSMIC signatures to mutational data.

    Arguments:
      -i, --input INPUT                      Path to mutational matrix (SBS/DBS/ID) file or folder of files
      -o, --output OUTPUT                    Output directory. Defaullt: output
      -s, --signature_type                   Type of signatures (SBS,DBS,ID). Default: SBS
      -c, --context_type                     Matrix format (e.g. SBS96, SBS288, DBS78, ID83, ID415). Default: SBS96/DBS78/ID83
      -g, --genome_type                      Choose from exome or genome data
      -d, --signature_database               Optional path to .txt file to include only selected signatures
      -e, --exclude_signature_subgroups      List of signature subgroups you don't want to analyze

    '''
    # create parser
    parser = argparse.ArgumentParser(
        description='Assign COSMIC signatures to mutation data')
    
    # Add parser arguments
    parser.add_argument('-i','--input', required=True, help='Path to mutational matrix (SBS/DBS/ID) file or folder of files')
    parser.add_argument('-o','--output', help='Output directory. Default: output', default='output')
    parser.add_argument('-s','--signature_type', choices=['SBS', 'DBS', 'ID'], default='SBS', help='Type of signatures. Default: SBS')
    parser.add_argument('-c','--context_type', help='Matrix format (e.g. SBS96, SBS288, DBS78, ID83, ID415). Default: SBS96/DBS78/ID83')
    parser.add_argument('-g','--genome_type', choices=['exome', 'genome'], default='genome', help='Choose from exome or genome data. Default: genome')
    parser.add_argument('-d','--signature_database', help='Optional path to .txt file to include only selected signatures', default=None)
    parser.add_argument('-e','--exclude_signature_subgroups', help='List of signature subgroups you don\'t want to analyze', default=None)
    
    args = parser.parse_args()
    gen_ex = False if args.genome_type == 'genome' else True

    #checking format
    default_context = {'SBS': 'SBS96', 'DBS': 'DBS78', 'ID': 'ID83'}

    valid_context = {
        'SBS': ['SBS96', 'SBS288', 'SBS1536', 'SBS6144'],
        'DBS': ['DBS78', 'DBS1248'],
        'ID':  ['ID28', 'ID83', 'ID415', 'ID862']
    }

    if args.context_type:
        if args.context_type.isdigit():
            context_type = f"{args.signature_type}{args.context_type}"
        else:
            context_type = args.context_type
        if context_type not in valid_context[args.signature_type]:
            raise ValueError(
                f"The selected format {context_type} is not compatible with the signature type {args.signature_type}. "
                f"Allowed formats: {valid_context}")
    else:
        context_type = default_context[args.signature_type]

    # checking file(s) and its (their) format
    input_path = args.input
    if os.path.isfile(input_path):
        files_to_process = [input_path] if os.path.splitext(input_path)[1].lower() == '.all' else []
    elif os.path.isdir(input_path):
        files_to_process = [
            os.path.join(input_path, f)
            for f in os.listdir(input_path)
            if os.path.isfile(os.path.join(input_path, f)) and os.path.splitext(os.path.join(input_path, f))[1].lower() == '.all'
        ]
    else:
        raise ValueError(f"Provided {input_path} is not a file or directory")
    
    if not files_to_process:
        raise ValueError(f"Provided {input_path} does not contain any suitable file(s)")
    else:
        print(f"Found {len(files_to_process)} file(s) to process: {', '.join([os.path.basename(f) for f in files_to_process])}")
        context_type = args.signature_type
        tasks = []
        for file_path in files_to_process:
            output_dir = os.path.join(args.output, os.path.splitext(os.path.basename(file_path))[0])
            tasks.append((file_path, output_dir, context_type, args.signature_database, gen_ex, args.exclude_signature_subgroups))

        # Multiprocessing analysis
        with Pool(processes=cpu_count()) as pool:
            list(tqdm(pool.imap_unordered(analyze, tasks), total=len(tasks), desc="Analysing"))
        print(f"Analysis completed. Results saved in {args.output}")

if __name__ == "__main__":
    main()
