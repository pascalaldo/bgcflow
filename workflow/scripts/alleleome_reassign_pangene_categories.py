import argparse
from pathlib import Path

import pandas as pd


def Gene_class_def(freq, threshold_99, threshold_15):
    """
    Classify gene based on frequency.

    Parameters:
    freq (int): Frequency of the gene.
    threshold_99 (int): Threshold for core genes.
    threshold_15 (int): Threshold for rare genes.

    Returns:
    str: Gene class ('Core', 'Rare', or 'Accessory').
    """
    if freq >= threshold_99:
        return "Core"
    elif freq < threshold_15:
        return "Rare"
    else:
        return "Accessory"


def process_roary_data(gp_binary_path, summary_path, output_file):
    """
    Process Roary data.

    Parameters:
    gp_binary_path (str): Binary gene presence file.
    summary_path (str): Original pangenome summary file.
    output_file (str): File to store output.
    """
    output_file = Path(output_file)
    # output_file.parent.mkdir(parents=True, exist_ok=True)

    # Loading the presence/absence matrix
    PAM = pd.read_csv(gp_binary_path, index_col=0)

    summary = pd.read_csv(summary_path)
    n_genome = len(PAM.columns)
    threshold_99 = round(n_genome * 0.99, 0)
    threshold_15 = round(n_genome * 0.15, 0)
    summary["pangenome_class_2"] = summary["No. isolates"].apply(
        Gene_class_def, threshold_99=threshold_99, threshold_15=threshold_15
    )
    summary.to_csv(output_file, index=False)


def main():
    """
    Main function to parse arguments and call the processing function.
    """
    parser = argparse.ArgumentParser(description="Process Roary data.")
    parser.add_argument(
        "--gp_binary", type=str, required=True, help="File containing the binary gene precence matrix."
    )
    parser.add_argument(
        "--summary", type=str, required=True, help="File containing the original pangenome summary file."
    )
    parser.add_argument("--output-file", required=True, help="File to store output")
    args = parser.parse_args()

    process_roary_data(args.gp_binary, args.summary, args.output_file)


if __name__ == "__main__":
    main()
