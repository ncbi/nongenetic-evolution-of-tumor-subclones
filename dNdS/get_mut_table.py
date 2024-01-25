import trisicell as tsc
import pandas as pd

##########
### Create mutation tables required for dndscv
##########

regime = "all"

# Sublines in each regime
if regime == "har":
    sublines =  ["C18", "C15", "C11", "C16"]
elif regime == "has":
    sublines = ["C1", "C4", "C22"]
elif regime == "las":
    sublines = ["C3", "C10", "C14"]
elif regime == "all":
    sublines = ["C1", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11",
                "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20",
                "C21", "C22", "C23", "C24"]
else: 
    print("Invalid regime.")
    exit()

# Load data
bwes = tsc.datasets.sublines_bwes()
df = bwes.var

# Take only SNV data
df = df[(df["kind"] == "synonymous SNV") | (df["kind"] == "nonsynonymous SNV")]

# Split off ensemble id suffix
df["ensemble_id"] = df.index.str.split(".").str[0]

# Get mutations
tsc_output = bwes.to_df(layer="trisicell_output").transpose()
tsc_output = tsc_output[~(tsc_output == 1).all(axis=1)]
tsc_output["ensemble_id"] = tsc_output.index.str.split(".").str[0]

# Take only sublines of interest and where they are not all zero
tsc_subset = tsc_output[tsc_output.columns[tsc_output.columns.isin(sublines + ["ensemble_id"])]]
tsc_subset = tsc_subset[~(tsc_subset[sublines] == 0).all(axis=1)]

# Take the SNV data for the ids we have mutation data for
df = df[df["ensemble_id"].isin(tsc_subset["ensemble_id"])]


# Take only the columns of interest, reformat, and save
df = df[["chrom", "position", "reference", "alteration"]]
df["chrom"] = df["chrom"].str[3:]
df["sampleID"] = "sample1"
df = df.rename(columns={"chrom": "chr", "position": "pos", "reference": "ref", "alteration":"mut"})
df.to_csv("dNdS/mut_tables/{}_snv_only.csv".format(regime))