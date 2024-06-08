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

def ecdna_freq_genes(aa_df, id_field, gene_field, oncogene_field, copy_num_field,
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
        oncogenes_only (optional bool, default=False): if the user wants to only
            see oncogene frequencies, this can be set to True
            
    Returns:
        Dictionary with key='GeneID', value=['Gene Frequency', set(ccles)].
        Access frequency of gene X by genes[X][0]
        Access ccles of gene X by genes[X][1]
    """    
    # Extract all rows with ecDNA classification into a new df
    ecdna = aa_df.loc['ecDNA']
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

def group_cell_lines(meta_df, achilles_field, group_field):
    """
    Groups cell-lines by a certain field in the metadata.

    Args:
        meta_df (DataFrame): metadata DataFrame generated from read_metadata()
        achilles_field (str): name of field in the metadata with Achilles IDs
        group_field (str): name of field in the metadata with the group of interest's classification

    Returns:
        A dictionary of the groups as keys and the values as Achilles IDs
    """
    group_dict = meta_df.groupby(group_field)[achilles_field].apply(list).to_dict()
    return group_dict

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
        if (gene not in esen_df):
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

        plt.savefig("data/Outputs/Gene Plots/" + nametag + "_" + gene + "_gene_plot.png")
        plt.close()

    print("Number of ecDNA+ cell lines:", len(ec_plus))
    print("Number of ecDNA- cell lines:", len(ec_minus))

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by="p-value")
    results_df.to_csv("data/Outputs/Gene Tables/" + nametag + "_gene_table.csv", index=False)
    return results_df

def cell_line_plot(gene_dict, esen_df, all_ids, ecDNA_present_ids, nametag):
    """
    Given a dict of genes (where values are "UP", "DOWN", or NaN regulation) and a list of cell-lines,
    generates a plot for each cell-line such that essentiality is along the y-axis for each of the 
    three regulation categories as columns.  
    
    Args:
        gene_dict ([str]): dict of genes as keys and regulation directions as values
        esen_df (DataFrame): matrix of CRISPR essentiality data
        all_ids (set()): set containing all Achilles IDs
        ecDNA_present_ids (set()): set containing all Achilles IDs with 
        nametag (str): nametag to differentiate various runs; added to output files
            
    Returns:
        DataFrame for each cell-line's statistics of the Kruskal-Wallis H test
    """    
    results = []

    for id in all_ids:
        # Skipping if ID not found in CRISPR data
        if (id not in esen_df.index):
            continue

        if (id in ecDNA_present_ids):
            presence = "present"
        else:
            presence = "absent"

        up_scores = []
        down_scores = []
        same_scores = []

        for gene in gene_dict.keys():
            # Skipping if gene not found in CRISPR data
            if (gene not in esen_df):
                continue

            esen_val = esen_df[gene][id]
            # Checking for NaN input, zz
            if esen_val == "nan" or esen_val == None or esen_val == "" or esen_val == "NaN" or pd.isna(esen_val):
                continue

            regulation = gene_dict[gene]

            if regulation == "UP":
                up_scores.append(esen_val)
            elif regulation == "DOWN":
                down_scores.append(esen_val)
            else:
                same_scores.append(esen_val)
            
        stat, p_val = scipy.stats.kruskal(up_scores, down_scores, same_scores)
        results.append({'Cell-Line': id, 'H-stat': stat, 'p-value': p_val})

        # Create a DataFrame from the input arrays
        data = {
            'Cell_Type': ['Upregulated'] * len(up_scores) + ['Downregulated'] * len(down_scores) + ['No Regulation (Same)'] * len(same_scores),
            'ecDNA_Score': up_scores + down_scores + same_scores
        }
        df = pd.DataFrame(data)

        # Create the box and whisker plot with adjusted width and neutral fill color
        plt.figure(figsize=(8, 10))
        sns.boxplot(x='Cell_Type', y='ecDNA_Score', data=df, showfliers=False, width=0.4,
                    boxprops=dict(facecolor="white", edgecolor="gray"),
                    whiskerprops=dict(color="gray"), capprops=dict(color="gray"), medianprops=dict(color="gray"))

        # Plot the data points
        sns.stripplot(x='Cell_Type', y='ecDNA_Score', data=df, jitter=True, color='black', edgecolor='gray', alpha=0.1)
        plt.title('Cell line = ' + id + "; ecDNA " + presence +  " (" + nametag + " genes)")
        plt.xlabel('ecDNA Regulation')
        plt.ylabel('CRISPR Essentiality Score')
        plt.tight_layout()
        plt.savefig("data/Outputs/Cell Line Plots/" + nametag + "_" + id + "_ecDNA_" + presence + "_cell_line_plot.png")
        plt.close()
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by="p-value")

    # Adding ecDNA presence information for each cell line's name
    def modify_cell_line(cell_line):
        if cell_line in ecDNA_present_ids:
            return f"{cell_line} (ecDNA present)"
        else:
            return f"{cell_line} (ecDNA absent)"
    results_df["Cell-Line"] = results_df["Cell-Line"].apply(modify_cell_line)
    results_df.to_csv("data/Outputs/Cell Line Tables/" + nametag + "_gene_table.csv", index=False)
    return results_df

def grouped_cell_line_plot(gene_dict, esen_df, group_ids, nametag):
    """
    Given a dict of genes (where values are "UP", "DOWN", or NaN regulation) and a list of grouped cell-lines,
    generates a plot for each category (group) such that essentiality is along the y-axis for each of the 
    three regulation columns.  
    
    Args:
        gene_dict ([str]): dict of genes as keys and regulation directions as values
        esen_df (DataFrame): matrix of CRISPR essentiality data
        group_ids (set()): set containing group IDs
        nametag (str): nametag to differentiate various runs; added to output files
            
    Returns:
        DataFrame for each cell-line's statistics of the Kruskal-Wallis H test
    """    
    results = []

    for id in group_ids:
        # Skipping if ID not found in CRISPR data
        if (id not in esen_df.index):
            continue

        up_scores = []
        down_scores = []
        same_scores = []

        for gene in gene_dict.keys():
            # Skipping if gene not found in CRISPR data
            if (gene not in esen_df):
                continue

            esen_val = esen_df[gene][id]
            # Checking for NaN input, zz
            if esen_val == "nan" or esen_val == None or esen_val == "" or esen_val == "NaN" or pd.isna(esen_val):
                continue

            regulation = gene_dict[gene]

            if regulation == "UP":
                up_scores.append(esen_val)
            elif regulation == "DOWN":
                down_scores.append(esen_val)
            else:
                same_scores.append(esen_val)
            
        stat, p_val = scipy.stats.kruskal(up_scores, down_scores, same_scores)
        results.append({'Cell-Line': id, 'H-stat': stat, 'p-value': p_val})

        # Create a DataFrame from the input arrays
        data = {
            'Cell_Type': ['Upregulated'] * len(up_scores) + ['Downregulated'] * len(down_scores) + ['No Regulation (Same)'] * len(same_scores),
            'ecDNA_Score': up_scores + down_scores + same_scores
        }
        df = pd.DataFrame(data)

        # Create the box and whisker plot with adjusted width and neutral fill color
        plt.figure(figsize=(8, 10))
        sns.boxplot(x='Cell_Type', y='ecDNA_Score', data=df, showfliers=False, width=0.4,
                    boxprops=dict(facecolor="white", edgecolor="gray"),
                    whiskerprops=dict(color="gray"), capprops=dict(color="gray"), medianprops=dict(color="gray"))

        # Plot the data points
        sns.stripplot(x='Cell_Type', y='ecDNA_Score', data=df, jitter=True, color='black', edgecolor='gray', alpha=0.1)
        plt.title('Cell line = ' + id + " (" + nametag + " genes)")
        plt.xlabel('ecDNA Regulation')
        plt.ylabel('CRISPR Essentiality Score')
        plt.tight_layout()
        plt.savefig("data/Outputs/Oncolineage Plots/" + nametag + "_" + id + "_oncolineage_cell_line_plot.png")
        plt.close()
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by="p-value")
    results_df.to_csv("data/Outputs/Oncolineage Tables/" + nametag + "_gene_table.csv", index=False)
    return results_df