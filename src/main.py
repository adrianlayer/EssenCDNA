import sys, argparse
from utils import execute_workflow

"""
    File:          main.py
    Description:   Main script to parse arguments & execute the EssenCDNA workflow

    Names:         Adrian Layer, Anthony Vasquez, Yasmin A. Jaber, Omar Halawa
    Emails:        alayer@ucsd.edu, pavasquez@ucsd.edu, yjaber@ucsd.edu, ohalawa@ucsd.edu
    Project:       EssenCDNA (final project for BENG/CSE/BIMM 182 @ UC San Diego)
    Repository:    https://github.com/adrianlayer/EssenCDNA
"""

# Function to parse command-line inputs
def parse_args():
    parser = argparse.ArgumentParser()

    # TODO: add parameters as necessary

    # File inputs
    parser.add_argument("--esen", help="DepMap CRISPR essentiality input", type=str, required=True)
    parser.add_argument("--aa", help="Amplicon Architect aggregated results data input", type=str, required=True)
    parser.add_argument("--meta", help="DepMap metadata input", type=str, required=True)

    # File field names
    parser.add_argument("--aa_type_field", help="AA type field in AA aggregated results data", type=str, required=True, defualt="Classification")
    parser.add_argument("--dep_ccle_field", help="CCLE name field in DepMap metadata", type=str, required=True, default="CCLEName")
    parser.add_argument("--aa_samples_field", help="AA cell-line field in AA aggregated results data", type=str, required=True, default="Sample name")
    parser.add_argument("--aa_gene_field", help="all genes field in AA aggregated results data", type=str, required=True, default="All genes")
    parser.add_argument("--aa_oncogene_field", help="oncogenes field in AA aggregated results data", type=str, required=True, default="Oncogenes")

    # Developer arguments
    parser.add_argument("-v", "--verbose", help="Verbosity flag", action='store_true')
    parser.add_argument("-d", "--debug", help="Debugging flag", action='store_true')

    args = parser.parse_args()
    return args

# Main function to execute workflow using parsed arguments
def main():
    args = parse_args()
    execute_workflow(args)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
        sys.exit(1)