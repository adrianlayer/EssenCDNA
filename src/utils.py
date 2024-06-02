import pandas as pd
from stats import *

"""
    File:          utils.py
    Description:   Script containing helper functions for io, plots, and executing workflow

    Names:         Adrian Layer, Anthony Vasquez, Yasmin A. Jaber, Omar Halawa
    Emails:        alayer@ucsd.edu, pavasquez@ucsd.edu, yjaber@ucsd.edu, ohalawa@ucsd.edu
    Project:       EssenCDNA (final project for BENG/CSE/BIMM 182 @ UC San Diego)
    Repository:    https://github.com/adrianlayer/EssenCDNA
"""

def read_esen(esen_filename):
    """
    Reads in DepMap CRISPR essentiality (dependency) data and outputs a standardized DataFrame of it.

    Args:
        esen_filename (str): filename of the DepMap CRISPR essentiality data (cell-lines by genes)

    Returns:
        DataFrame of the DepMap essentiality data
    """
    crispr_df = pd.read_csv(esen_filename, index_col=0)

    # Trimming off gene names
    crispr_df.columns = [
        col.split()[0] for col in crispr_df.columns
    ]

    return crispr_df

def read_metadata(metadata_filename):
    """
    Reads in DepMap metadata and outputs a DataFrame of it.

    Args:
        metadata_filename (str): filename of the metadata

    Returns:
        DataFrame of the DepMap metadata
    """
    meta_df = pd.read_csv(metadata_filename, index_col=0)
    return meta_df

def ccle_to_achilles(meta_df, ccle_field):
    """
    Creates a dictionary for cell-line name conversion from CCLE name to the Achilles ID.

    Args:
        meta_df (DataFrame): metadata DataFrame generated from read_metadata()
        ccle_field (str): name of field in the metadata with CCLE names

    Returns:
        Dictionary for name (ID) conversion
    """

    id_convert = {y: x for x, y in zip(meta_df.index.tolist(), meta_df[ccle_field].tolist())}
    return id_convert

def execute_workflow(args):
    return

# print(read_esen("../data/CRISPRGeneEffect.csv"))
# print(read_esen("../data/CRISPRGeneEffect.csv")["A1BG"]["ACH-000001"])
# print(read_metadata("../data/Model.csv"))
# print(ccle_to_achilles(read_metadata("../data/Model.csv"), "CCLEName")["59M_OVARY"])