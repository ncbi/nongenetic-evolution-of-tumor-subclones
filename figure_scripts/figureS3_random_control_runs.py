import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns 

df = pd.read_csv("results/control/results.csv")

sns.histplot(df["total"])

ax = plt.gca()
ax.vlines(616, ymin=0, ymax=1, transform=ax.get_xaxis_transform(), color="red")
ax.vlines(812, ymin=0, ymax=1, transform=ax.get_xaxis_transform(), color="red")
ax.vlines(1277, ymin=0, ymax=1, transform=ax.get_xaxis_transform(), color="red")

ymax = ax.get_ylim()[1]
plt.text(x=616+10, y=ymax-1, s="MA-S: 616")
plt.text(x=812+10, y=ymax-1, s="HA-R: 812")
plt.text(x=1277-140, y=ymax-1, s="SA-S: 1,277")

plt.title("Distribution of the number of genes predicted to have\nadaptive expression when 3 random sublines are tested")
plt.xlabel("# genes predicted to have adaptive expression")
plt.ylabel("# of random runs")
plt.tight_layout()
plt.savefig("figures/cdist.svg")