import pandas as pd
from collections import defaultdict
from ast import literal_eval
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import scipy
from scipy import stats

"""
    File:          utils.py
    Description:   Script containing helper functions for io, plots, and executing workflow

    Names:         Adrian Layer, Anthony Vasquez, Yasmin A. Jaber, Omar Halawa
    Emails:        alayer@ucsd.edu, pavasquez@ucsd.edu, yjaber@ucsd.edu, ohalawa@ucsd.edu
    Project:       EssenCDNA (final project for BENG/CSE/BIMM 182 @ UC San Diego)
    Repository:    https://github.com/adrianlayer/EssenCDNA
"""

###############################################################################
# Functions to parse input files

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

def read_aa(aa_results_filename, gene_id_field):
    """
    Reads in AA results and outputs a DataFrame of it.

    Args:
        aa_results_filename (str): filename of the AA aggregated results data
        aa_type_field (str): field name containing amplicon type

    Returns:
        DataFrame of the AA results
    """
    # Indexed by amplicon classification
    aa_df = pd.read_csv(aa_results_filename, index_col=gene_id_field)
    return aa_df

def read_literature_gene_list(list_filename, gene_field, direction_field):
    """
    Reads in gene list from literature ("ecDNA Target genes.csv) and outputs a dictionary of it
        (key is gene, value is ecDNA regulation direction).

    Args:
        list_filename (str): filename of ecDNA list with gene vs regulation direction
        gene_field (str): field name for gene names (Gene|Gene ID)
        direction_field (str): field name for regulation direction

    Returns:
        Dictionary of gene (keys) and regulation (values)
    """

    list_df = pd.read_csv(list_filename)

    def extract_gene_name(value):
        if isinstance(value, str):
            return value.split('|')[0]
        return 'NA'

    list_df['Gene'] = list_df[gene_field].apply(extract_gene_name)

    # Create the dictionary
    gene_dict = pd.Series(list_df[direction_field].values, index=list_df['Gene']).to_dict()
    del gene_dict["NA"]
    return gene_dict

###############################################################################
# Functions to manipulate data

def avg_esen(gene_list, crispr_df):
    """
    Reads in a list of genes and a CRISPR essentiality DataFrame and outputs a DataFrame of average essentiality values.

    Args:
        gene_list ([str]]): list of genes as strings
        crispr_df (DataFrame): matrix of CRISPR essentiality values

    Returns:
        DataFrame with each gene and its mean, median, and stdev essentiality values
    """
    data = []

    for gene in gene_list:
        if gene in crispr_df:
            mean_esen = crispr_df[gene].mean()
            median_esen = crispr_df[gene].median()
            stdev_esen = crispr_df[gene].std()

            if median_esen > 0:
                category = "Non-essential"
            elif median_esen < -1:
                category = "Common-essential"
            else:
                category = "Selectively-essential"

            data.append([gene, mean_esen, median_esen, stdev_esen, category])

    # Creating a DataFrame from the data
    columns = ['Gene', 'Mean', 'Median', 'Standard Deviation', 'Category']
    result_df = pd.DataFrame(data, columns=columns)
    return result_df

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

def group_achilles(meta_df, ):
    """
    Groups Achilles IDs (cell-lines) based on primary disease, onco-lineage, and onco-subtype.

    Args:
        meta_df (DataFrame): metadata DataFrame generated from read_metadata()
        ccle_field (str): name of field in the metadata with CCLE names

    Returns:
        Three dictionaries, each with keys representing categories and values as Achilles IDs. 
    """
    oncotree_primary_disease_dict = meta_df.groupby('OncotreePrimaryDisease')['ModelID'].apply(list).to_dict()
    oncotree_lineage_dict = meta_df.groupby('OncotreeLineage')['ModelID'].apply(list).to_dict()
    oncotree_subtype_dict = meta_df.groupby('OncotreeSubtype')['ModelID'].apply(list).to_dict()

    return oncotree_primary_disease_dict, oncotree_lineage_dict, oncotree_subtype_dict

def ecdna_freq_genes(aa_df, id_field, gene_field, oncogene_field, copy_num_field, classification_field,
                     oncogenes_only: bool=False) -> dict[str, list[float, set[str]]]:
    """
    Creates a descending sorted dictionary of all genes or oncogenes by frequency
    across all ccles. Genes are mapped to a length 2 list of [frequency, 
    {set of all ccles gene is a part of}]
    
    Args:
        aa_df (DataFrame): Amplicon Architect data matrix
        id_field (str): field name for samples names (or IDs); typically is first column
        gene_field (str): field name for all genes in an amplicon
        oncogene_field (str): field name for all oncogenes in an amplicon
        copy_num_field (str): field name for median copy count of feature in an amplicon
        classification_field (str): field name for classification of ccle feature in an amplicon
        oncogenes_only (optional bool, default=False): if the user wants to only
            see oncogene frequencies, this can be set to True
            
    Returns:
        Dictionary with key='GeneID', value=['Gene Frequency', set(ccles)].
        Access frequency of gene X by genes[X][0]
        Access ccles of gene X by genes[X][1]
    """
    ec_aa_df = aa_df.set_index(classification_field)
    
    # Extract all rows with ecDNA classification into a new df
    ecdna = ec_aa_df.loc['ecDNA']
    # Apply literal eval so that the gene_field column does not get evaluated as str
    ecdna.loc[:,oncogene_field] = ecdna[oncogene_field].apply(literal_eval) # use .loc to avoid chain indexing
    ecdna.loc[:,gene_field] = ecdna[gene_field].apply(literal_eval)

    def def_val():
        """
        Helper for frequent genes dictionary definition.
        """
        return[0.0, set()]

    # Use lambda function to set default gene count to 0
    genes = defaultdict(def_val)

    for i, row in ecdna.iterrows():
        if oncogenes_only:
            for gene in row[oncogene_field]:
                genes[gene][0] += row[copy_num_field]
                genes[gene][1].add(row[id_field])
        else:
            for gene in row[gene_field]:
                genes[gene][0] += row[copy_num_field]
                genes[gene][1].add(row[id_field])
    sorted_genes = {k: v for k, v in sorted(genes.items(), key=lambda x: x[1][0], reverse=True)}
    
    return sorted_genes

###############################################################################
# Functions for plotting graphs and generating statistics

def plot_genes(gene_list, esen_df, ec_aa_plus_ids, nametag):
    """
    Given a list of genes and a list of cell-lines that are ecDNA+, generates a plot for each gene 
    such that essentiality is along the y-axis for each of the two categories as columns (ecDNA +/-).  
    
    Args:
        gene_list ([str]): list of genes as strings
        esen_df (DataFrame): matrix of CRISPR essentiality data
        ec_aa_plus_ids (set()): set containing Achilles IDs that contain ecDNA
        nametag (str): nametag to differentiate various runs; added to output files
            
    Returns:
        DataFrame for each gene's statistics of the Mann-Whitney U rank test
    """
    results = []

    for gene in gene_list:
        
        # Checking for possibility of not finding the gene in CRISPR knockout data
        if gene not in esen_df:
            continue

        ec_plus = []
        ec_minus = []

        for index, row in esen_df.iterrows():
            if index in ec_aa_plus_ids:
                ec_plus.append(esen_df[gene][index])
            else:
                ec_minus.append(esen_df[gene][index])

        # # Calculate median values
        # median_ec_plus = np.median(ec_plus)
        # median_ec_minus = np.median(ec_minus)

        stat, p_val = scipy.stats.mannwhitneyu(ec_plus, ec_minus)
        results.append({'Gene': gene, 'u-stat': stat, 'p-value': p_val})

        # Create a DataFrame from the input arrays
        data = {
            'Cell_Type': ['ecDNA+'] * len(ec_plus) + ['ecDNA-'] * len(ec_minus),
            'ecDNA_Score': ec_plus + ec_minus
        }

        df = pd.DataFrame(data)

        # Create the box and whisker plot with adjusted width
        plt.figure(figsize=(8, 10))
        sns.boxplot(x='Cell_Type', y='ecDNA_Score', data=df, showfliers=False, width=0.4, showmeans=False)
        sns.stripplot(x='Cell_Type', y='ecDNA_Score', data=df, color='black', alpha=0.1, size=4, jitter=True, edgecolor='none')
        plt.title('Gene = ' + gene)
        plt.xlabel('ecDNA Presence')
        plt.ylabel('CRISPR Essentiality Score')
        plt.tight_layout()

        # # Add annotations for median values
        # plt.text(0, median_ec_plus, f'Median: {median_ec_plus:.2f}', horizontalalignment='center', verticalalignment='bottom')
        # plt.text(1, median_ec_minus, f'Median: {median_ec_minus:.2f}', horizontalalignment='center', verticalalignment='bottom')

        plt.savefig("../data/Outputs/Gene Plots/" + nametag + "_" + gene + "_gene_plot.png")
        plt.close()

    print("Number of ecDNA+ cell lines:", len(ec_plus))
    print("Number of ecDNA- cell lines:", len(ec_minus))

    results_df = pd.DataFrame(results)
    results_df.to_csv("../data/Outputs/Gene Tables/" + nametag + "_gene_table.csv", index=False)
    return results_df

def execute_workflow(args):
    """
    Executes overall pipeline through calling helper functions.
    
    Args:
        args (argparse.Namespace): parsed arguments from command-line

    Returns:
        None
    """    
    esen = read_esen(args.esen)
    meta = read_metadata(args.meta)
    aa = read_aa(args.aa, args.aa_type_field)
    ecdna_list = read_ecdna_list(args)

    id_dict = ccle_to_achilles(meta, args.dep_ccle_field)


    return None

###############################################################################
# Quick debugging code
# print(read_esen("../data/CRISPRGeneEffect.csv"))
# print(read_esen("../data/CRISPRGeneEffect.csv")["A1BG"]["ACH-000001"])
# print(read_metadata("../data/Model.csv"))
# print(ccle_to_achilles(read_metadata("../data/Model.csv"), "CCLEName")["59M_OVARY"])
# print(ecdna_freq_genes(read_aa("../data/aggregated_results.csv","Classification"), "Sample name", "All genes", "Oncogenes", True))
# print((read_ecdna_list("../data/ecDNA Target genes.csv", "#Gene|GeneId", "direction(for_geneset_enrichment)")))
# print(avg_esen(list(read_ecdna_list("../data/ecDNA Target genes.csv", "#Gene|GeneId", "direction(for_geneset_enrichment)").keys()), read_esen("../data/CRISPRGeneEffect.csv")))