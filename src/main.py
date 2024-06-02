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