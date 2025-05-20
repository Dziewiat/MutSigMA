import os
import pandas as pd
from tqdm import tqdm

def filter_database(project_names: list[str] | None = None,
                    patient_ids: list[str] | None = None,
                    primary_sites: list[str] | None = None,
                    chromosomes: list[str] | None = None,
                    database_filepath: str = 'data/mutations.parquet.gzip') -> pd.DataFrame:
    '''Filter mutations data based on user requests (list format).'''
    # Load mutations database
    print('Loading mutational database...')
    data = pd.read_parquet(database_filepath)
    
    # Filter database by requested parameters
    print('Locating requested records...')
    if project_names:
        data = data.loc[data['project_short_name'].isin(project_names)]
    if patient_ids:
        data = data.loc[data['case_barcode'].isin(patient_ids)]
    if primary_sites:
        data = data.loc[data['primary_site'].isin(primary_sites)]
    if chromosomes:
        data = data.loc[data['Chromosome'].isin(chromosomes)]

    # Check if there are any records matching all the parameters
    assert not data.empty, 'There are no records matching the chosen parameters. Check for typing errors in the request.'

    return data


def extract_vcf(data: pd.DataFrame):
    '''Extract and save VCF files from dataframe into a specified folder.'''
    # Convert dataframe into vcf-friendly format
    df_vcf = data.apply(
        lambda row: [
            row['Chromosome'].replace('chr', ''),
            row['Start_Position'],
            row['case_barcode'],
            row['Reference_Allele'],
            row['Tumor_Seq_Allele2'],
            '.',
            'Simulations',
            'GRCh38',
            row['Reference_Allele'],
            '+1'
        ],
        axis=1
    )

    # Create new VCF dataframe
    columns_vcf = [
        "Chromosome", "Position", "ID", "Ref", "Alt", "Qual", "Filter", "Info", "Motif", "Strand"
    ]
    df_vcf = pd.DataFrame(df_vcf.tolist(), columns=columns_vcf)

    # Group by 'case_barcode' (patient ID)
    grouped = df_vcf.groupby("ID")

    # Save each patient's data into a separate VCF file
    output_dir = "data/VCF"  # directory for VCF storage
    os.makedirs(output_dir, exist_ok=True)

    vcf_file_paths = []

    for file_gdc_id, group in tqdm(grouped, desc='VCF files extraction...'):
        vcf_file_name = f"{file_gdc_id}.vcf"
        vcf_file_path = os.path.join(output_dir, vcf_file_name)
        save_to_vcf(group, vcf_file_path)
        vcf_file_paths.append(vcf_file_path)


def save_to_vcf(group, file_path):
    '''Save patient data into a VCF file.'''

    # Create a VCF file header
    vcf_header = """##fileformat=VCFv4.2
##source=CustomConversionScript
##reference=GRCh38
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tMOTIF\tSTRAND
"""

    # Write every patient's VCF file
    with open(file_path, "w") as vcf_file:
        vcf_file.write(vcf_header)
        for _, row in group.iterrows():
            vcf_file.write(
                f"{row['Chromosome']}\t{row['Position']}\t{row['ID']}\t"
                f"{row['Ref']}\t{row['Alt']}\t{row['Qual']}\t{row['Filter']}\t{row['Info']}\t"
                f"{row['Motif']}\t{row['Strand']}\n"
            )
