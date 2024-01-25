import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm
import seaborn as sns

##########
### Creates a draft of the cluster heatmap seen in figure 1c
### Modifications were made in Adobe Illustrator
##########

cluster_file = "results/tpm/kmeans_clusters.csv"
save_file = "figures/figure1c_cluster_heatmap_code_output.svg"
regimes = ["har", "has", "las"]

df = pd.read_csv(cluster_file)
df = df.set_index("gene_name", drop=True)
df = df[regimes + ["clusters"]]  # this will order regimes as we want
df = df.sort_values("clusters")

clusts = pd.DataFrame(df["clusters"])
df = df.drop(["clusters"], axis=1)
cmap = cm.get_cmap("tab20", 12)

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True, gridspec_kw={'width_ratios': [3, 1]})
sns.heatmap(data = df, ax=ax1, cmap="seismic",  cbar_kws = {"location":"left"})
sns.heatmap(data = clusts, ax=ax2, cmap=cmap) # adding this so we can define cluster boundaries
plt.yticks([])
plt.savefig(save_file)
plt.show()