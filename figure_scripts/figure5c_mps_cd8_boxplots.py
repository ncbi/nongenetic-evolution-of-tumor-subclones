import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import stats

save_path = "figures/treatment_data/mps_cd8.svg"

data = pd.read_csv("data/mouse_treatment/treatment_mps_cd8.csv")
data = data.replace("Stable", "Responder") # Include stable disease in responders
print(data)

mps_t = stats.ttest_ind(data[data["Response"] == "Responder"]["Adjusted MPS"], data[data["Response"] == "No-response"]["Adjusted MPS"])
cd8_t = stats.ttest_ind(data[data["Response"] == "Responder"]["CD8"], data[data["Response"] == "No-response"]["CD8"])
print("MPS t-test:", mps_t)
print("CD8 t-test:", cd8_t)

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
sns.boxplot(data=data, x="Response", y="Adjusted MPS", ax=ax1)
ax1.set_title("Adjusted MPS")
sns.boxplot(data=data, x="Response", y="CD8", ax=ax2)
ax2.set_title("CD8 score")
plt.tight_layout()
plt.savefig(save_path)
plt.show()