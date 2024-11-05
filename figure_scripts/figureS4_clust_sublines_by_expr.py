from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
from scipy.stats import spearmanr

# Cluster sublines based on expression

data_path = "data/tpm/"
result_path = "results/tpm/gene_lists/adpt_{}/all_results.csv"
regime_path = "regime_files/"

har_sublines = ["C11", "C15", "C16", "C18"]
sas_sublines = ["C1", "C4", "C22"]
mas_sublines = ["C3", "C10", "C14"]

def color_labels(ax):
    labels = ax.get_xmajorticklabels()
    for label in labels:
        if label.get_text() in har_sublines:
            label.set_color("blue")
        if label.get_text() in sas_sublines:
            label.set_color("orange")
        if label.get_text() in mas_sublines:
            label.set_color("brown")

print("Loading data.")
# Expression data
data_df = pd.DataFrame()
for f in os.listdir(data_path):
    new_df = pd.read_csv(data_path+f)
    new_df["gene"] = new_df["gene"].str.split("_", expand=True)[1]
    data_df = pd.concat([data_df, new_df])
data_avg_df = data_df.drop("replicate", axis=1)
data_avg_df = data_avg_df.groupby(["gene", "species"]).mean().reset_index()
data_avg_df = data_avg_df.pivot(columns="gene", index="species", values="exprval")

# Adaptivity data
har_df = pd.read_csv(result_path.format("har"))[["gene_name", "qvalue"]]
sas_df = pd.read_csv(result_path.format("sas"))[["gene_name", "qvalue"]]
mas_df = pd.read_csv(result_path.format("mas"))[["gene_name", "qvalue"]]

har_adpt_5 = har_df[har_df["qvalue"] < 0.05]["gene_name"].to_numpy()
sas_adpt_5 = sas_df[sas_df["qvalue"] < 0.05]["gene_name"].to_numpy()
mas_adpt_5 = mas_df[mas_df["qvalue"] < 0.05]["gene_name"].to_numpy()

har_adpt_1 = har_df[har_df["qvalue"] < 0.01]["gene_name"].to_numpy()
sas_adpt_1 = sas_df[sas_df["qvalue"] < 0.01]["gene_name"].to_numpy()
mas_adpt_1 = mas_df[mas_df["qvalue"] < 0.01]["gene_name"].to_numpy()

adpt_5_np = np.union1d(np.union1d(har_adpt_5, sas_adpt_5), mas_adpt_5)
adpt_1_np = np.union1d(np.union1d(har_adpt_1, sas_adpt_1), mas_adpt_1)

# Non-adaptive is always > 0.05, regardless of what the p-value cutoff is for adaptive
har_na = har_df[har_df["qvalue"] > 0.05]["gene_name"].to_numpy()
sas_na = sas_df[sas_df["qvalue"] > 0.05]["gene_name"].to_numpy()
mas_na = mas_df[mas_df["qvalue"] > 0.05]["gene_name"].to_numpy()
na_np = np.union1d(np.union1d(har_na, sas_na), mas_na)

method="ward"
metric="euclidean"

fig, axes = plt.subplots(nrows=1, ncols=4, sharey=True, figsize=(15, 6))

# 1. Tree with all
temp = data_avg_df.to_numpy()
labels = data_avg_df.index.to_list()
Z = linkage(temp, method=method, metric=metric)
Z[:,2] = np.arange(5, 115, 5) # Make the branches even
dendrogram(Z, labels=labels, color_threshold=0, ax=axes[0])
color_labels(axes[0])
axes[0].set_title("All genes")

# 2. Tree wtih neutral genes
temp = data_avg_df[na_np].to_numpy()
labels = data_avg_df.index.to_list()
Z = linkage(temp, method=method, metric=metric)
Z[:,2] = np.arange(5, 115, 5) # Make the branches even
dendrogram(Z, labels=labels, color_threshold=0, ax=axes[1])
color_labels(axes[1])
axes[1].set_title("Non-adaptive genes")

# 3. Tree with adaptive genes p < 0.05
temp = data_avg_df[adpt_5_np].to_numpy()
labels = data_avg_df.index.to_list()
Z = linkage(temp, method=method, metric=metric)
Z[:,2] = np.arange(5, 115, 5) # Make the branches even
R = dendrogram(Z, labels=labels, color_threshold=0, ax=axes[2])
color_labels(axes[2])
axes[2].set_title("Genes with adaptive expression q<{}".format(0.05))

# 4. Tree with adaptive genes p < 0.01
temp = data_avg_df[adpt_1_np].to_numpy()
labels = data_avg_df.index.to_list()
Z = linkage(temp, method=method, metric=metric)
Z[:,2] = np.arange(5, 115, 5) # Make the branches even
R = dendrogram(Z, labels=labels, color_threshold=0, ax=axes[3])
color_labels(axes[3])
axes[3].set_title("Genes with adaptive expression q<{}".format(0.01))

save_file = "figures/figureS4_subline_clusts.svg"
fig.suptitle("Ward clustering of sublines by expression")
plt.tight_layout()
plt.savefig(save_file)
plt.show()