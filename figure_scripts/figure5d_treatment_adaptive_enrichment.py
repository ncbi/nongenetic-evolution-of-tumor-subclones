import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

df_save_file = "results/mouse_treatment/differential_expression/deseq2_summary.csv"
plt_save_file = "figures/figure3d_adaptive_in_diff_eq_code_output.svg"

# DEGs that are evaluated for adaptivity
deg_adpt = pd.read_csv("results/mouse_treatment/differential_expression/deseq2_adpt_merge.csv")
deg_adpt["r_nr"] = ["R > NR" if g["log2FoldChange"] > 0 else "R < NR" for i, g in deg_adpt.iterrows()]

r = len(deg_adpt[deg_adpt["log2FoldChange"] > 0])
nr = len(deg_adpt[deg_adpt["log2FoldChange"] < 0])

df_data = []

for regime in ["har", "has", "las"]:
    curr_deg = deg_adpt[deg_adpt["{} q-value".format(regime)] < 0.05]
    
    r_pos = len(curr_deg[(curr_deg["log2FoldChange"] > 0) & (curr_deg["{} logFC".format(regime)] > 0)])
    r_neg = len(curr_deg[(curr_deg["log2FoldChange"] > 0) & (curr_deg["{} logFC".format(regime)] < 0)])
    nr_pos = len(curr_deg[(curr_deg["log2FoldChange"] < 0) & (curr_deg["{} logFC".format(regime)] > 0)])
    nr_neg = len(curr_deg[(curr_deg["log2FoldChange"] < 0) & (curr_deg["{} logFC".format(regime)] < 0)])

    df_data += [[regime, "R > NR", "positive", r_pos, r_pos/r]]
    df_data += [[regime, "R > NR", "negative", r_neg, r_neg/r]]
    df_data += [[regime, "R < NR", "positive", nr_pos, nr_pos/nr]]
    df_data += [[regime, "R < NR", "negative", nr_neg, nr_neg/nr]]

df = pd.DataFrame(df_data, columns = ["regime","r_nr","pos_neg_adpt","adpt_de","adpt_de/de"])

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, sharey=True)
df["adpt_de/de"] = round(df["adpt_de/de"] * 100, 1)
print(df)
df.to_csv(df_save_file)

# Plot, separating R > N and R < N
b1 = sns.barplot(data=df[df["r_nr"] == "R > NR"], x="regime", y="adpt_de/de", hue="pos_neg_adpt", ax=ax1)
b2 = sns.barplot(data=df[df["r_nr"] == "R < NR"], x="regime", y="adpt_de/de", hue="pos_neg_adpt", ax=ax2)
ax1.set_title("R > N")
ax2.set_title("R < N")

# Create a single legend
handles = []
labels = []
axHandle, axLabel = ax1.get_legend_handles_labels()
handles.extend(axHandle)
labels.extend(axLabel)

# Create single labels
ax1.set_xlabel("")
ax1.set_ylabel("")
ax2.set_ylabel("")
fig.supylabel("Enrichment of adaptive genes")

# Add title
fig.suptitle("Enrichment of Adaptive Genes\nin Genes Differentially Expressed\nBetween Responders and Non-Responders\n(DESeq2)")

plt.tight_layout()
plt.savefig(plt_save_file)
plt.show()