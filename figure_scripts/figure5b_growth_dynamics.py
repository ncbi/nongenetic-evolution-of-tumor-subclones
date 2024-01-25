import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

data = pd.read_csv("data/mouse_treatment/post_treatment_growth_dynamics.csv")
data = data.rename(columns = {"Unnamed: 0": "day"})

data = data.melt(id_vars = "day", var_name="mouse")
data["response"] = data["mouse"].str.split(".", expand=True)[0]
print(data)

colors = {"anti-CTLA-4 (Non-responder)": "blue", "anti-CTLA-4 (Responder)": "orange", "IgG2b": "gray"}

f, (ax1, ax2) = plt.subplots(2, 1, figsize=(5,7), sharex=True, gridspec_kw={'height_ratios': [1, 3]})

for mouse in data["mouse"].unique():
    curr_data = data[data["mouse"] == mouse]["response"].reset_index(drop=True)
    color = colors[curr_data[0]]
    sns.lineplot(data = data[data["mouse"] == mouse], x="day", y="value", color=color, ax=ax1)
    sns.lineplot(data = data[data["mouse"] == mouse], x="day", y="value", color=color, ax=ax2)

# Split axis
ax1.set_ylim(400, 1000)  # outliers only
ax2.set_ylim(0, 300)  # most of the data
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
d = .015  # how big to make the diagonal lines in axes coordinates
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal


plt.savefig("figures/treatment_data/growth_dynamics_code_output.svg")
plt.show()