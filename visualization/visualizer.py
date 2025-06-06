import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
import os
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")


SIGNATURE_ETIOLOGY = {
    # Aging
    'SBS1': 'Clock-like', 'SBS5': 'Clock-like',
    # Smoking
    'SBS4': 'Tobacco', 'SBS29': 'Tobacco', 'SBS92': 'Tobacco', 'DBS2': 'Tobacco', 'ID3': 'Tobacco',
    # UV exposure
    'SBS7a': 'UV', 'SBS7b': 'UV', 'SBS7c': 'UV', 'SBS7d': 'UV',
    'SBS38': 'UV', 'DBS1': 'UV', 'ID13': 'UV',
    # DNA repair deficiency
    'SBS3': 'HRD', 'DBS13': 'HRD', 'ID6': 'HRD',
    'SBS6': 'dMMR', 'SBS14': 'dMMR', 'SBS15': 'dMMR', 'SBS20': 'dMMR', 'SBS21': 'dMMR',
    'SBS26': 'dMMR', 'SBS44': 'dMMR', 'DBS7': 'dMMR', 'DBS10': 'dMMR', 'ID7': 'dMMR',
    'SBS30': 'dBER',
    'ID8': 'NHEJ',
    # Polymerase
    'SBS9': 'Polymerase mutation',  'SBS10a': 'Polymerase mutation', 'SBS10b': 'Polymerase mutation',
    'SBS10c': 'Polymerase mutation',  'SBS10d': 'Polymerase mutation', 'DBS3': 'Polymerase mutation',
    # Slippage
    'ID1': 'Replication slippage', 'ID2': 'Replication slippage',
    # AID/APOBEC
    'SBS2': 'AID/APOBEC activity', 'SBS13': 'AID/APOBEC activity', 'SBS84': 'AID/APOBEC activity',
    'SBS85': 'AID/APOBEC activity',
    # Chemotherapy
    'SBS11': 'Chemotherapy', 'SBS25': 'Chemotherapy', 'SBS31': 'Chemotherapy', 'SBS32': 'Chemotherapy',
    'SBS35': 'Chemotherapy', 'SBS86': 'Chemotherapy', 'SBS87': 'Chemotherapy', 'SBS90': 'Chemotherapy',
    'SBS99': 'Chemotherapy', 'DBS5': 'Chemotherapy',
    # Aristolochic acid
    'SBS22a': 'Aristolochic acid', 'SBS22b': 'Aristolochic acid', 'DBS20': 'Aristolochic acid', 'ID23': 'Aristolochic acid',
    # Aflatoxin
    'SBS24': 'Aflatoxin',
    # Colibactin
    'SBS88': 'Colibactin', 'ID18': 'Colibactin',
    # Haloalkanes
    'SBS42': 'Haloalkanes',
    # ROS
    'SBS17b': 'ROS', 'SBS18': 'ROS', 'SBS36': 'ROS',
    # Topoisomerase activity
    'ID4': 'TOP mutation', 'ID17': 'TOP mutation',
    # Unknown
    'SBS8': 'Unknown', 'SBS12': 'Unknown', 'SBS16': 'Unknown', 'SBS17a': 'Unknown',
    'SBS19': 'Unknown', 'SBS23': 'Unknown', 'SBS27': 'Unknown', 'SBS28': 'Unknown',
    'SBS33': 'Unknown', 'SBS34': 'Unknown', 'SBS37': 'Unknown', 'SBS39': 'Unknown',
    'SBS40a': 'Unknown', 'SBS40b': 'Unknown', 'SBS40c': 'Unknown', 'SBS41': 'Unknown',
    'SBS43': 'Unknown', 'SBS45': 'Unknown', 'SBS46': 'Unknown', 'SBS47': 'Unknown',
    'SBS48': 'Unknown', 'SBS49': 'Unknown', 'SBS50': 'Unknown', 'SBS51': 'Unknown',
    'SBS52': 'Unknown', 'SBS53': 'Unknown', 'SBS54': 'Unknown', 'SBS55': 'Unknown',
    'SBS56': 'Unknown', 'SBS57': 'Unknown', 'SBS58': 'Unknown', 'SBS59': 'Unknown',
    'SBS60': 'Unknown', 'SBS89': 'Unknown', 'SBS91': 'Unknown', 'SBS93': 'Unknown',
    'SBS94': 'Unknown', 'SBS95': 'Unknown', 'SBS96': 'Unknown', 'SBS97': 'Unknown',
    'SBS98': 'Unknown', 'DBS4': 'Unknown', 'DBS6': 'Unknown', 'DBS8': 'Unknown',
    'DBS9': 'Unknown', 'DBS11': 'Unknown', 'DBS12': 'Unknown', 'DBS14': 'Unknown',
    'DBS15': 'Unknown', 'DBS16': 'Unknown', 'DBS17': 'Unknown', 'DBS18': 'Unknown',
    'DBS19': 'Unknown', 'ID5': 'Unknown', 'ID9': 'Unknown', 'ID10': 'Unknown',
    'ID11': 'Unknown', 'ID12': 'Unknown', 'ID14': 'Unknown', 'ID15': 'Unknown',
    'ID16': 'Unknown', 'ID19': 'Unknown', 'ID20': 'Unknown', 'ID21': 'Unknown',
    'ID22': 'Unknown',
}

def find_input_file(input_path):
    """ Find the Assignment_Solution_Activities.txt file """
    if os.path.isfile(input_path):
        return input_path, extract_dataset_name(input_path)

    activities_file = os.path.join(input_path, "Assignment_Solution", "Activities",
                                   "Assignment_Solution_Activities.txt")

    if os.path.exists(activities_file):
        dataset_name = extract_dataset_name(input_path)
        return activities_file, dataset_name
    else:
        print(f"Error: Could not find Assignment_Solution_Activities.txt in {input_path}")
        print(f"Expected path: {activities_file}")
        return None, None

def extract_dataset_name(path):
    """ Extract dataset name from path (e.g. MutSigMA.SBS96 from output/MutSigMA.SBS96/) """
    path_parts = Path(path).parts
    return path_parts[-1]

def load_data(input_file):
    """ Load mutational signatures data """
    try:
        data = pd.read_csv(input_file, sep='\t', index_col=0)
        data = data.apply(pd.to_numeric, errors='coerce').fillna(0)
        print(f"Successfully loaded data: {data.shape[0]} samples, {data.shape[1]} signatures")
        return data
    except Exception as e:
        print(f"Error loading data: {e}")
        return None

def calculate_figure_width(n_signatures, min_width=12, max_width=50):
    """ Calculate optimal figure width based on number of signatures """
    width = max(min_width, min(max_width, n_signatures * 0.5))
    return width

def create_boxplot(data, output_dir, dataset_name, show_outliers=True, show_only=False):
    """ Create boxplot of signature activities """
    # Filter out signatures with all zeros
    active_data = data.loc[:, (data != 0).any(axis=0)]
    n_signatures = active_data.shape[1]

    if n_signatures == 0:
        print("No active signatures found for boxplot.")
        return

    width = calculate_figure_width(n_signatures)
    fig, ax = plt.subplots(figsize=(width, 7))

    data_melted = active_data.melt(var_name='Signature', value_name='Activity')

    if show_outliers:
        sns.boxplot(data=data_melted, x='Signature', y='Activity', ax=ax)
        ax.set_title('Mutational Signatures Activity (with outliers)', fontweight='bold', fontsize=16)
        filename = f"{dataset_name}_boxplot_with_outliers.png"
    else:
        sns.boxplot(data=data_melted, x='Signature', y='Activity', showfliers=False, ax=ax)
        ax.set_title('Mutational Signatures Activity (no outliers)', fontweight='bold', fontsize=16)
        filename = f"{dataset_name}_boxplot_no_outliers.png"

    ax.set_xlabel('Mutational Signatures', labelpad=20, fontsize=14)
    ax.set_ylabel('Activity Level', labelpad=20, fontsize=14)
    ax.tick_params(axis='x', rotation=90)

    if show_only:
        plt.subplots_adjust(bottom=0.2)
        plt.show()
    else:
        plt.tight_layout()
        output_path = os.path.join(output_dir, filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Boxplot saved: {output_path}")
    plt.close()


def create_barplot_active_signatures(data, output_dir, dataset_name, show_only=False):
    """ Create barplot showing number of patients with active signatures """
    # Count patients with active signatures (activity > 0)
    active_counts = (data > 0).sum()
    active_counts = active_counts[active_counts > 0]

    if len(active_counts) == 0:
        print("No active signatures found for barplot.")
        return

    width = calculate_figure_width(len(active_counts))
    fig, ax = plt.subplots(figsize=(width, 7))

    bars = ax.bar(range(len(active_counts)), active_counts.values)
    ax.set_title('Number of Patients with Active Signatures', fontweight='bold', fontsize=16)
    ax.set_xlabel('Mutational Signatures', labelpad=20, fontsize=14)
    ax.set_ylabel('Number of Patients', labelpad=20, fontsize=14)
    ax.set_xticks(range(len(active_counts)))
    ax.set_xticklabels(active_counts.index, rotation=90)

    # Value labels
    for i, v in enumerate(active_counts.values):
        ax.text(i, v + 0.2, str(v), ha='center', va='bottom')

    if show_only:
        plt.subplots_adjust(bottom=0.2)
        plt.show()
    else:
        plt.tight_layout()
        output_path = os.path.join(output_dir, f"{dataset_name}_active_signatures_barplot.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Active signatures barplot saved: {output_path}")
    plt.close()

def create_etiology_piechart(data, output_dir, dataset_name, sample_id=None, show_only=False):
    """ Create pie chart of signature etiologies based on prevalence or signature count """
    if sample_id and sample_id not in data.index:
        print(f"Error: Sample {sample_id} not found in data.")
        return

    etiology_counts = {}

    if sample_id:
        # For single patient: count active signatures per etiology
        patient_data = data.loc[sample_id]
        active_signatures = patient_data[patient_data > 0].index

        for sig in active_signatures:
            etiology = SIGNATURE_ETIOLOGY.get(sig, 'Unknown')
            etiology_counts[etiology] = etiology_counts.get(etiology, 0) + 1

        title_text = f'Active Signature Etiologies - Patient {sample_id}'
        filename = f"{dataset_name}_etiology_piechart_{sample_id}.png"
        unit_text = "signatures"

    else:
        # For all patients: count in how many patients each etiology appears
        for etiology_name in set(SIGNATURE_ETIOLOGY.values()) | {'Unknown'}:
            etiology_signatures = [sig for sig, etiol in SIGNATURE_ETIOLOGY.items() if etiol == etiology_name]
            if etiology_name == 'Unknown':
                all_mapped_sigs = set(SIGNATURE_ETIOLOGY.keys())
                unknown_sigs = set(data.columns) - all_mapped_sigs
                etiology_signatures.extend(unknown_sigs)

            # Count patients with at least one active signature with this etiology
            etiology_signatures_in_data = [sig for sig in etiology_signatures if sig in data.columns]
            if etiology_signatures_in_data:
                patients_with_etiology = (data[etiology_signatures_in_data] > 0).any(axis=1).sum()
                if patients_with_etiology > 0:
                    etiology_counts[etiology_name] = patients_with_etiology

        title_text = 'Etiology Prevalence - All Patients'
        filename = f"{dataset_name}_etiology_piechart_overall.png"
        unit_text = "patients"

    if not etiology_counts:
        print("No active signatures found for pie chart.")
        return

    # Sorted legend
    sorted_etiologies = sorted(etiology_counts.items(), key=lambda x: x[1], reverse=True)
    sorted_values = [item[1] for item in sorted_etiologies]

    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_position([0.1, 0.1, 0.5, 0.8])

    wedges, texts = ax.pie(sorted_values, startangle=90)

    total = sum(sorted_values)
    if sample_id:
        legend_labels = [f"{etiology} ({count} {unit_text}, {count / total * 100:.1f}%)"
                         for etiology, count in sorted_etiologies]
    else:
        legend_labels = [f"{etiology} ({count} {unit_text}, {count / total * 100:.1f}%)"
                         for etiology, count in sorted_etiologies]

    ax.legend(wedges, legend_labels, title="Etiologies",
              loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

    ax.set_title(title_text)

    if show_only:
        plt.show()
    else:
        plt.tight_layout()
        output_path = os.path.join(output_dir, filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Etiology pie chart saved: {output_path}")
    plt.close()

def cluster_signatures(data, output_dir=None, dataset_name=None, method='complete', show_only=False):
    """Cluster signatures hierarchically and visualize as a heatmap with dendrograms."""
    
    # Filter out signatures with all zeros
    active_data = data.loc[:, (data != 0).any(axis=0)]
    if active_data.shape[1] == 0:
        print("No active signatures found for clustering heatmap.")
        return
    
    # Remove .txt extension if present
    clean_name = dataset_name.replace('.txt', '') if dataset_name else 'heatmap'

    # Create clustermap (z-score normalization by signature for better visualization)
    cg = sns.clustermap(
        active_data,
        metric='cosine',
        method=method,
        cmap='viridis',
        figsize=(max(12, active_data.shape[1] * 0.5), 10),
        z_score=1,  # normalize by signature (column)
        linewidths=0.1,
        cbar_kws={'label': 'Z-score Activity'},
        cbar_pos=(1.02, 0.2, 0.01, 0.6)
    )

    cg.ax_heatmap.set_xlabel('Mutational Signatures', fontsize=14)
    cg.ax_heatmap.set_ylabel('Samples', fontsize=14)
    cg.ax_heatmap.set_title('Hierarchical Clustering of Mutational Signatures', fontsize=16, fontweight='bold')

    if show_only:
        plt.subplots_adjust(right=0.90)
        plt.show()
    else:
        if output_dir is None or dataset_name is None:
            print("Output directory and dataset name required to save clustering heatmap.")
            return
        output_path = os.path.join(output_dir, f"{clean_name}_signature_clustermap.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Signature clustering heatmap saved: {output_path}")
    plt.close()

def generate_summary_report(data, output_dir, dataset_name):
    """ Generate a summary report of the analysis """
    report_path = os.path.join(output_dir, f"{dataset_name}_summary_report.txt")

    with open(report_path, 'w') as f:
        f.write("MutSigMA - Mutational Signatures Analysis Summary Report\n")
        f.write("=" * 57 + "\n\n")

        f.write(f"Dataset: {dataset_name}\n")
        f.write(f"Dataset Overview:\n")
        f.write(f"- Number of samples: {data.shape[0]}\n")
        f.write(f"- Number of signatures: {data.shape[1]}\n\n")

        # Number of patients with each signature
        prevalence = (data > 0).sum().sort_values(ascending=False)
        f.write("Top 10 Most Prevalent Signatures (by number of patients):\n")
        for i, (sig, count) in enumerate(prevalence.head(10).items(), 1):
            etiology = SIGNATURE_ETIOLOGY.get(sig, 'Unknown')
            percentage = (count / data.shape[0]) * 100
            f.write(f"{i:2d}. {sig}: {count} patients ({percentage:.1f}%) - {etiology}\n")
        f.write("\n")

        # Etiology distribution
        etiology_counts = {}
        for sig, activity in data.sum().items():
            if activity > 0:
                etiology = SIGNATURE_ETIOLOGY.get(sig, 'Unknown')
                etiology_counts[etiology] = etiology_counts.get(etiology, 0) + 1

        f.write("Etiology Distribution:\n")
        f.write("-" * 40 + "\n")
        f.write("By number of different signatures:\n")
        for etiology, count in sorted(etiology_counts.items(), key=lambda x: x[1], reverse=True):
            f.write(f"- {etiology}: {count} signatures\n")

        f.write("\nBy patient prevalence (how many patients have each etiology):\n")
        etiology_prevalence = {}
        for etiology_name in set(SIGNATURE_ETIOLOGY.values()) | {'Unknown'}:
            etiology_signatures = [sig for sig, etiol in SIGNATURE_ETIOLOGY.items() if etiol == etiology_name]
            if etiology_name == 'Unknown':
                all_mapped_sigs = set(SIGNATURE_ETIOLOGY.keys())
                unknown_sigs = set(data.columns) - all_mapped_sigs
                etiology_signatures.extend(unknown_sigs)

            # Count patients with at least one active signature from this etiology
            etiology_signatures_in_data = [sig for sig in etiology_signatures if sig in data.columns]
            if etiology_signatures_in_data:
                patients_with_etiology = (data[etiology_signatures_in_data] > 0).any(axis=1).sum()
                if patients_with_etiology > 0:
                    etiology_prevalence[etiology_name] = patients_with_etiology

        for etiology, count in sorted(etiology_prevalence.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / data.shape[0]) * 100
            f.write(f"- {etiology}: {count} patients ({percentage:.1f}%)\n")

    print(f"Summary report saved: {report_path}")


def main():
    '''
    MutSigMA Visualization Tool

    Arguments:
      -i, --input                            Path to Assignment_Solution folder
      -o, --output                           Output directory. Default: plots
      -b, --boxplot                          Generate boxplot of signature activities
      -n, --no_outliers                      Hide outliers in boxplot
      -r, --barplot                          Generate barplot of active signatures count
      -p, --piechart                         Generate pie chart of signature etiologies
      -d, --id                               Specific sample ID for individual pie chart
      -f, --first_n                          Generate pie charts for first N patients (default: 0)
      -a, --all                              Generate all visualizations
      -t, --report                           Generate summary report
      -s, --show                             Display plots only (no saving)
      -c, --cluster_signatures               Generate clustering heatmap of signatures
    '''

    parser = argparse.ArgumentParser(description=
                                     'MutSigMA Visualization Tool')

    parser.add_argument('-i', '--input', required=True, help='Path to Assignment_Solution folder')
    parser.add_argument('-o', '--output', default='plots', help='Output directory (default: plots)')
    parser.add_argument('-b', '--boxplot', action='store_true', help='Generate boxplot of signature activities')
    parser.add_argument('-n', '--no_outliers', action='store_true', help='Hide outliers in boxplot')
    parser.add_argument('-r', '--barplot', action='store_true', help='Generate barplot of active signatures count')
    parser.add_argument('-p', '--piechart', action='store_true', help='Generate pie chart of signature etiologies')
    parser.add_argument('-d', '--id', help='Specific sample ID for individual pie chart')
    parser.add_argument('-f', '--first_n', type=int, default=0, help='Generate pie charts for first N patients (default: 0)')
    parser.add_argument('-a', '--all', action='store_true', help='Generate all visualizations')
    parser.add_argument('-t', '--report', action='store_true', help='Generate summary report')
    parser.add_argument('-s', '--show', action='store_true', help='Display plots only (no saving)')
    parser.add_argument('-c', '--cluster_signatures', action='store_true', help='Generate clustering heatmap of signatures')

    args = parser.parse_args()

    # Find input file
    input_file, dataset_name = find_input_file(args.input)
    if input_file is None:
        sys.exit(1)

    print(f"Using input file: {input_file}")
    print(f"Dataset name: {dataset_name}")

    # Load data
    data = load_data(input_file)
    if data is None:
        sys.exit(1)

    # Create output directory if not --show arg
    if not args.show:
        os.makedirs(args.output, exist_ok=True)
        print(f"Output directory: {args.output}")

    if args.all or args.boxplot:
        print("\nGenerating boxplot...")
        if args.all:
            # If -all, generate both
            create_boxplot(data, args.output, dataset_name, show_outliers=True, show_only=args.show)
            create_boxplot(data, args.output, dataset_name, show_outliers=False, show_only=args.show)
        else:
            # If --boxplot: only with outliers, if --boxplot --no_outliers: only no outliers
            create_boxplot(data, args.output, dataset_name, show_outliers=not args.no_outliers, show_only=args.show)

    if args.all or args.barplot:
        print("\nGenerating barplot...")
        create_barplot_active_signatures(data, args.output, dataset_name, show_only=args.show)

    if args.all or args.piechart:
        print("\nGenerating pie chart...")

        # If --piechart --first_n: generate only individual patient charts
        if args.first_n > 0:
            sample_list = data.index[:args.first_n]
            for sample in sample_list:
                create_etiology_piechart(data, args.output, dataset_name, sample_id=sample, show_only=args.show)
        elif args.id:
            # If --piechart --id: generate piechart for specific sample ID
            create_etiology_piechart(data, args.output, dataset_name, sample_id=args.id, show_only=args.show)
        else:
            # If --piechart: generate chart for all patients (default)
            create_etiology_piechart(data, args.output, dataset_name, sample_id=None, show_only=args.show)
    if args.all or args.cluster_signatures:
        print("\nGenerating cluster heatmap...")
        cluster_signatures(data, args.output, dataset_name, show_only=args.show)
        
    if (args.all or args.report) and not args.show:
        print("\nGenerating summary report...")
        generate_summary_report(data, args.output, dataset_name)

    if not any([args.boxplot, args.barplot, args.piechart, args.all, args.report, args.cluster_signatures]):
        print("No visualization option selected. Use --help for available options.")
        print("Quick start: python visualizer.py -i output_folder/output_file --all")


if __name__ == "__main__":
    main()