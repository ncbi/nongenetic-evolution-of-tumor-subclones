import pandas as pd
import os
from scipy import stats
from statsmodels.stats.multitest import multipletests
import numpy as np

##########
### Run differential expression over the simulated data
### Update file paths to reflect your data and save locations
##########

base_path = "data/simulated/neut_sim/"
regime ="har"
regime_file = "regime_files/{}.csv".format(regime)
save_path = "results/simulated/differential_expression/neut_sim/"

def run_de(df, r1_taxa, r2_taxa, save_file):
    df["chosen"] = df["species"].isin(r2_taxa)

    r1_data = df[df["species"].isin(r1_taxa)]
    r2_data = df[df["species"].isin(r2_taxa)]

    results = []

    for gene in df["gene"].unique():
        ensemble_id = gene
        gene_name = gene
        curr_data = df[df["gene"] == gene]
        r1_data = curr_data[curr_data["species"].isin(r1_taxa)]
        r2_data = curr_data[curr_data["species"].isin(r2_taxa)]
        
        r1_mean = r1_data["exprval"].mean()
        r2_mean = r2_data["exprval"].mean()
        r1_std = r1_data["exprval"].std()
        r2_std = r2_data["exprval"].std()

        ind_t = stats.ttest_ind(r1_data["exprval"].to_list(), r2_data["exprval"].to_list())
        logFC = np.log2(r2_mean/r1_mean)
        results.append([ensemble_id, gene_name, r1_mean, r2_mean, r2_mean-r1_mean, r1_std, r2_std, logFC, ind_t[1]])

    results_df = pd.DataFrame(data=results, columns=["ensemble_id", "gene_name", "base_mean", "adpt_mean", "mean_diff(adpt-base)", "base_std", "adpt_std", "logFC", "tt_pvalue"])
    fdr_p=multipletests(pvals=results_df["tt_pvalue"], alpha=0.05, method="fdr_bh")
    results_df["fdr_p"] = fdr_p[1]

    print("# diff expr:", len(results_df[results_df["fdr_p"] < 0.05]))
    if save_path:
        results_df.to_csv(save_file, index=False)


regime_df = pd.read_csv(regime_file)

r1_taxa = regime_df[(regime_df["regime"] == "background") &  (regime_df["node2"].isna())]["node"].to_list()
r2_taxa = regime_df[(regime_df["regime"] == "chosen") &  (regime_df["node2"].isna())]["node"].to_list()

df = pd.DataFrame()

for f in os.listdir(base_path):
    df = pd.read_csv(base_path + f)
    if save_path:
        run_de(df, r1_taxa, r2_taxa, save_path + f[:-4] +  ".csv")
    else:
        run_de(df, r1_taxa, r2_taxa, save_path)
