import trisicell as tsc
import pandas as pd
import numpy as np
import math
import argparse

##########
## Save the Trisicell data into usable format.
## Keep only genes that have non-zero values for all species.
## Log normalize data.
##########

def main(start, end, inc, dtype, save_path, save_prefix):

    # Load the data
    sc_data = tsc.datasets.sublines_scrnaseq()
    exp_data = sc_data["expression"]

    while start < end:

        # Set boundary values for this gene set
        curr_end = start + inc
        if curr_end > end: curr_end = end

        # Get the gene indices for this file
        gene_indices = list(range(start,curr_end))
        gene_str = "{}-{}".format(start, curr_end)

        print("Running for genes " + gene_str)


        # Create new dataframe for this gene set
        exp_df = pd.DataFrame(columns=["gene", "species", "replicate", "exprval"])
        
        dropped = 0

        # Loop through the genes for this file
        for index in gene_indices:

            # Get the name of this gene (ensembleid_genename)
            gene_name = exp_data.var.index[index]

            # Loop through the sublines (C1, C2, etc)
            for subline in exp_data.obs["clone"].unique(): 
                
                # C2 isn't in our tree, so skip it
                if subline == "C2":
                    continue

                # Get the replicates (subclones) for this subline (C1_1, C1_2, etc)
                subclones = exp_data[exp_data.obs["clone"] == subline].obs.index.tolist()

                # Loop through the replicates (subclones)
                for subclone in subclones:

                    rep_name = "R{}".format(subclone.split("_")[-1])

                    # Get the expression data for all genes for this replicate
                    rep = exp_data[exp_data.obs.index == subclone].layers[dtype]

                    # Error check formatting
                    if len(rep) != 1:
                        print("error", subline)
                        exit()

                    # Get the raw expression value for this gene
                    rep = rep[0][index]

                    # If the value is 0, skip and move to the next
                    if rep == 0:
                        continue
                    else:
                        # If not zero, do log transformation
                        expr_val = math.log2(1+rep)

                    # Add this value to the dataframe
                    curr_df = pd.DataFrame.from_dict({
                        "gene": [gene_name],
                        "species": [subline],
                        "replicate": [rep_name],
                        "exprval": [expr_val]
                    })
                    exp_df = pd.concat([exp_df, curr_df], ignore_index = True)

            # If the gene has all NA values, remove it
            if (exp_df[exp_df["gene"] == gene_name]["exprval"] == "NA").all():
                # print("Gene {} is all NA - dropping.".format(gene_name))
                exp_df = exp_df.drop(exp_df[exp_df["gene"] == gene_name].index)
                dropped += 1

            # If the gene doesn't have a value for every species, remove it
            elif (len(exp_df[exp_df["gene"] == gene_name]["species"].unique()) < 23):
                exp_df = exp_df.drop(exp_df[exp_df["gene"] == gene_name].index)
                dropped += 1

        print("Dropped {} genes".format(dropped))

        # Save the current gene set
        exp_df.to_csv("{}{}gene_{}.csv".format(save_path, save_prefix, gene_str), index=False)

        # Increase indices to next gene set
        start = curr_end
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, action="store", default=0,
                        help="The index of the gene to start analysis with. Default = 0.")
    parser.add_argument("--end", type=int, action="store", default=55401,
                        help="The last gene to include in analysis. Default = 55401 (total number of genes in dataset).")
    parser.add_argument("--inc", type=int, action="store", default=1000,
                        help="The increment of genes to analyze and save separately. " +
                        "If genes are excluded because of non-zero values, the save files will " + 
                        "have fewer than inc genes. Default = 1000.")
    parser.add_argument("--dtype", type=str, action="store", default="tpm",
                        help="The normalization type of the data. Either tpm or fpkm. Default = tpm.")
    parser.add_argument("--save_path", type=str, action="store", default="",
                        help="The directory into which the data will be saved. Default = ./")
    parser.add_argument("--save_prefix", type=str, action="store", default="",
                        help="The prefix for the saved data files (for instance if saving multiple " +
                        "types of data to the save directory). Default = ''.")
    args = parser.parse_args()

    if not args.dtype in ["tpm", "fpkm"]:
        print("dtype must be either 'tpm' or 'fpkm'.")
        exit()
    if args.save_prefix != "" and args.save_prefix[-1] != "_":
        args.save_prefix += "_"

    main(args.start, args.end, args.inc, args.dtype, args.save_path, args.save_prefix)
