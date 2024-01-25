import trisicell as tsc
import pandas as pd
import numpy as np

##########
### Count the number of genes that are mutated and adaptive
### (Supplement Section 3)
##########

regime = "har"
# regime = "has"
# regime = "las"

base_path = "results/tpm/gene_lists/adpt_{}/".format(regime)
adpt_path = base_path + "all_results.csv"
save_path = "{}_{}.csv".format("mutation_info", regime)

print("Regime:", regime)

# Load mutation data
bwes = tsc.datasets.sublines_bwes()
tsc_output = bwes.to_df(layer="trisicell_output").transpose()
# Remove clonal mutations
tsc_output = tsc_output[~(tsc_output == 1).all(axis=1)]
# Remove genes without mutations
tsc_output = tsc_output[~(tsc_output == 0).all(axis=1)]
# Label by gene
tsc_output["ensemble_id"] = tsc_output.index.str.split(".").str[0]
# Group by gene
tsc_by_gene = tsc_output.groupby("ensemble_id").sum()
# Make mutations binary per gene rather than sum
tsc_by_gene[tsc_by_gene > 0] = 1
print("Total # of mutated genes:", len(tsc_by_gene))

# Load adaptive genes for selected regime
adpt_res_all = pd.read_csv(adpt_path)
adpt_res_all["ensemble_id"] = adpt_res_all["ensemble_id"].str.split(".").str[0]
adpt_res = adpt_res_all[adpt_res_all["adaptive"] == "adaptive"]
n_mut = tsc_by_gene[tsc_by_gene.index.isin(adpt_res_all["ensemble_id"])]
print("Total # of mutated genes evaluated for adaptivity:", len(n_mut))
print("Total # of adaptive genes:", len(adpt_res))

# Get genes that are adative and mutated
tsc_adpt = tsc_by_gene[tsc_by_gene.index.isin(adpt_res["ensemble_id"])]
print("# of genes both adaptive and mutated:", len(tsc_adpt))

print("# of genes that are mutated and adaptive / # genes that are mutated:", round(len(tsc_adpt) / len(n_mut), 3))