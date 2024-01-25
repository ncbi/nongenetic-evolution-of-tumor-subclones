import pandas as pd
import numpy as np
from sklearn import cluster
from sklearn.metrics import silhouette_score

##########
### Run k-means over the genes labeled adaptive in at least one regime
### Replace results_dir and save_file as needed
### Because of random seed, the clusters will be in the same order
### Genes may not be in the same order in the final data frame, but they will be
### ordered by cluster
##########

def get_num_clusters(df):
    sil = []
    kmin = 2
    kmax = 20

    for k in range(kmin, kmax):
        preds = cluster.KMeans(n_clusters=k, random_state=0).fit(df[["har", "has", "las"]])
        sil.append(silhouette_score(df[["har", "has", "las"]], preds.labels_, metric = "euclidean"))

    return np.argmax(sil) + 2 # add two because we start k at 2

# Only adaptive genes in these files
results_dir = "results/tpm/gene_lists/adpt_{}/gene_info.csv"
save_file = "results/tpm/kmeans_clusters.csv"
regimes = ["har", "has", "las"]

df = pd.DataFrame(columns=["gene_name", "regime", "v"])

for regime in regimes:
    temp = pd.read_csv(results_dir.format(regime))
    print(regime, len(temp))
    temp["regime"] = regime
    temp["v"] = np.where(temp["logFC"] < 0, np.log2(temp["qvalue"]), -1 * np.log2(temp["qvalue"])) # all q values will be < 1, so log will be negative, we want negative when log fold is negative
    df = pd.concat([df, pd.DataFrame(temp[["gene_name", "regime", "v"]])])
    # df = df.append(temp[["gene_name", "regime", "v"]])
# print(len(df))
# exit()

df = df.pivot(index = "gene_name", columns="regime", values="v")
df = df.fillna(0) # genes that aren't adaptive in that regime are assigned 0
print(len(df))

k = get_num_clusters(df)

preds = cluster.KMeans(n_clusters=k, random_state=0).fit_predict(df[["har", "has", "las"]])
df["clusters"] = preds

# We rename the clusters to sort then better for the heatmap -- cluster order was chosen manually
# Reordering will work on rerunning the script because of the set random state
x     = [0, 1, 2, 3, 4,  5, 6, 7, 8, 9, 10, 11]
order = [2, 5, 4, 7, 0, 11, 1, 8, 6, 10, 9, 3]
x = np.linspace(0, 11, 12)
order_pairs = list(zip(order, x))
for i in range(len(order_pairs)):
    df.loc[df["clusters"] == order_pairs[i][0],"order"] = order_pairs[i][1]
df["clusters"] = df["order"].astype(int)
df = df.drop("order", axis=1)

df = df.sort_values("clusters")
print(df)

df.to_csv(save_file)