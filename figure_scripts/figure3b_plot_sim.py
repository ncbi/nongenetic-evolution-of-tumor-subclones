from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats

##########
## Figure 3b
## Plots boxplots of the true positive and false positive predictons
## of EvoGeneX and differental expression analysis on simulated data
## across all parameter combinations
##########

sim_dist = "poisson" # poisson or nb
save_path = "figures/"
adpt_sim_adpt_run = pd.read_csv("results/{}_sim/adpt_sim_adpt_run.csv".format(sim_dist))
adpt_sim_de = pd.read_csv("results/{}_sim/adpt_sim_de.csv".format(sim_dist))
neut_sim_adpt_run = pd.read_csv("results/{}_sim/neut_sim_adpt_run.csv".format(sim_dist))
neut_sim_de = pd.read_csv("results/{}_sim/neut_sim_de.csv".format(sim_dist))

if sim_dist == "poisson":
    adpt_vars = ["sigma_sq", "alpha", "theta_ratio"]
    neut_vars = ["sigma_sq"]
else:
    adpt_vars = ["sigma_sq", "r", "alpha", "theta_ratio"]
    neut_vars = ["sigma_sq", "r"]    

adpt_sim_adpt_run.columns = adpt_vars + ["# predicted adaptive"]
neut_sim_adpt_run.columns = neut_vars + ["# predicted adaptive"]
adpt_sim_de.columns = adpt_vars + ["de"]
neut_sim_de.columns = neut_vars + ["de"]

res = stats.ttest_ind(adpt_sim_adpt_run["# predicted adaptive"], adpt_sim_de["de"])
print("adpt sim t-test:", res)
res = stats.ttest_ind(neut_sim_adpt_run["# predicted adaptive"], neut_sim_de["de"])
print("neut sim t-test:", res)

adpt_sim = adpt_sim_adpt_run.merge(adpt_sim_de, on=adpt_vars)
adpt_sim = pd.melt(adpt_sim, id_vars=adpt_vars, var_name="pred_type", value_name="count", ignore_index=True)
# print(adpt_sim)

neut_sim = neut_sim_adpt_run.merge(neut_sim_de, on=neut_vars)
neut_sim = pd.melt(neut_sim, id_vars=neut_vars, var_name="pred_type", value_name="count", ignore_index=True)
# print(neut_sim)

plt.figure()
sns.boxplot(data=adpt_sim, x="pred_type", y="count", hue="pred_type")
plt.title("Adaptively simulated data - true positives")
plt.savefig("{}figure3b_sim_true_positives_code_output.svg".format(save_path))
plt.show()
plt.close()

plt.figure()
sns.boxplot(data=neut_sim, x="pred_type", y="count", hue="pred_type")
plt.title("Neutrally simulated data - false positives")
plt.savefig("{}figure3b_sim_false_positives_code_output.svg")
plt.show()
plt.close()