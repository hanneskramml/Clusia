#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import io
from statsmodels.stats.multitest import multipletests
import argparse

def load_tsv(directory, file_ending, sep="\t", how="all"):
    files_list = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith(file_ending) == True:
                f = os.path.normpath(os.path.join(root, filename))
                files_list.append(f)
    dictionary = {}
    for file in files_list:
        key = file
        value = pd.read_csv(file, sep=sep).dropna(axis=0, how=how)
        dictionary[key] = value
    return dictionary

def load_tsv_list(directory, file_ending, sep="\t"):
    files_list = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith(file_ending) == True:
                f = os.path.normpath(os.path.join(root, filename))
                files_list.append(f)
    df_list = []
    for file in files_list:
        df = pd.read_csv(file, sep=sep)
        df_list.append(df)
    return df_list
                
def adjust_pvalues(p_values, method, alpha):
    p_value_bonf = multipletests(p_values, method=method, alpha=alpha)[1] 
    corr_alpha_bonf = multipletests(p_values, method=method, alpha=alpha)[3]
    return p_value_bonf, corr_alpha_bonf

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory", type=str, default="./fimo", help="Define root directory.")
parser.add_argument("-p", "--prefixes", nargs="+", help="Define species. Species names are folders \
                     containing fimo result files as well as parts of filenames. Species names \
                         are separated by a single space.")

args = parser.parse_args()

# load fimo results files
directory = args.directory
for prefix in args.prefixes:
    print(f"Working on {prefix} fimo reults ...")
    print("Loading input files...")
    fimos = load_tsv(f"{directory}/{prefix}", "fimo.tsv", sep="\t", how="any")
    
    # load transcript ids
    with open (f"{prefix}.cam_genes.txt", "r") as f:
        data = f.read()
        gene_ids = data.split("\n")
    
    # load dictionary with motif ID and function of the acossiated TF
    motifs_dict = {}
    with open(f"{directory}/TF_function.txt") as f:
        for line in f:
            (key, val) = line.split("\t")
            motifs_dict[key] = val.rstrip("\n") 
    
    # correct p-value using Bonferroni correction and filter significant results
    print("Correct p-values using Bonferroni correction...")
    try:
        for key, fimo in fimos.items():
            # fimo["description"] = fimo["description"].replace({"%2C": ","}, regex=True)
            fimo["sequence_name"] = fimo["sequence_name"].str.split("::", n=1).str[0]
            fimo = fimo.drop(["strand"], axis=1) 
            fimo.insert(loc=7, column="p_value_bonf", value=np.nan)
            fimo = fimo.rename({"p-value": "p_value", "q-value": "q_value",
                                "motif_id": "JASPAR_ID", "motif_alt_id": "motif_name", 
                                "sequence_name": "gene"}, axis=1)
            fimo["p_value_bonf"], corr_alpha_bonf = adjust_pvalues(fimo["p_value"], "bonferroni", 1)
            fimo = fimo.sort_values("p_value")
            # fimo["significance_after_correction"] = fimo["p_value"] <= corr_alpha_bonf
            fimos[key] = fimo
    except:
        pass
    
    # concatenate and fimo results
    fimos_concat = pd.concat([fimo for fimo in fimos.values()])
    bools = fimos_concat["p_value_bonf"] <= 0.05
    fimos_concat = fimos_concat[bools]
    fimos_concat.insert(loc=2, column="TF_function", value=np.nan)
    fimos_concat_copy = fimos_concat.copy()
    fimos_concat_copy.loc[:,"TF_function"] = fimos_concat["motif_name"].map(motifs_dict)
    fimos_concat = fimos_concat_copy
    fimos_concat = fimos_concat.sort_values("p_value")
    
    # get list of gene IDs for which no significant motifs were found
    no_motif = []
    for gene_id in gene_ids:
        try:  
            fimo_filtered = {key: fimo_run for key, fimo_run in fimos.items() if gene_id in key}
            fimo_df = pd.concat([fimo_run for key, fimo_run in fimo_filtered.items()])
            bools = fimo_df["p_value_bonf"] <= 0.05
            fimo_df = fimo_df[bools]
            if fimo_df.empty:
                # print(f"No significant motif occurences were found in {gene_id}.")
                no_motif.append(f"{gene_id}")
        except:
            no_motif.append(f"{gene_id}")
            # print(f"No motif occurences were found in {gene_id}.")
    
    # save list of non-signifincant motifs
    print("Saving list of gene IDs without motif occurences...")
    with io.open(f"{directory}/{prefix}.no_motif.txt", "w", newline="\n") as f:
        for element in no_motif:
            f.write(element + "\n")
            
    # save results as tsv
    print(f"Saving corrected fimo results as tsv file in {directory}...")
    fimos_concat.to_csv(f"{directory}/{prefix}.motifs.tsv", sep="\t", index=False)

# concatenate fimo result files of all species---------------------------------------------------------------
fimo_results = load_tsv_list(directory, ".motifs.tsv", sep="\t")
fimo_results_concat = pd.concat(fimo_results, axis=0, ignore_index=True)
fimo_results_concat.to_csv(f"{directory}/Clusia.motifs.tsv", sep="\t", index=False)

         