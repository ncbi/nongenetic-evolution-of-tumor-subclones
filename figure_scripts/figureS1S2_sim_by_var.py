from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats

##########
## Figures S1 and S2
## Plots boxplots of the true positive and false positive predictons
## of EvoGeneX and differental expression analysis on simulated data
## separating out the data by parameters
##########

sim_dist = "poisson" # poisson or nb

if sim_dist == "poisson":
    save_path = "figures/figureS1"
else:
    save_path = "figures/figureS2"

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

adpt_sim = adpt_sim_adpt_run.merge(adpt_sim_de, on=adpt_vars)
adpt_sim = pd.melt(adpt_sim, id_vars=adpt_vars, var_name="pred_type", value_name="count", ignore_index=True)
# print(adpt_sim)

neut_sim = neut_sim_adpt_run.merge(neut_sim_de, on=neut_vars)
neut_sim = pd.melt(neut_sim, id_vars=neut_vars, var_name="pred_type", value_name="count", ignore_index=True)
# print(neut_sim)

# Adaptive plots
for group_var in adpt_vars:
    print(group_var)

    # Test significance
    for val in adpt_sim[group_var].unique():
        res = stats.ttest_ind(adpt_sim[(adpt_sim[group_var] == val) & (adpt_sim["pred_type"] == "# predicted adaptive")]["count"],
                              adpt_sim[(adpt_sim[group_var] == val) & (adpt_sim["pred_type"] == "de")]["count"])
        print("{} t-test: {}".format(val, res))
    plt.figure()
    sns.boxplot(data=adpt_sim, x=group_var, y="count", hue="pred_type")
    plt.title(group_var)
    plt.savefig("{}_adpt_{}.svg".format(save_path, group_var))
    plt.show()
    plt.close()


# Neutral plots
for group_var in neut_vars:
    print(group_var)
    # Test significance
    for val in neut_sim[group_var].unique():
        res = stats.ttest_ind(neut_sim[(neut_sim[group_var] == val) & (neut_sim["pred_type"] == "# predicted adaptive")]["count"],
                              neut_sim[(neut_sim[group_var] == val) & (neut_sim["pred_type"] == "de")]["count"])
        print("{} t-test: {}".format(val, res))
    plt.figure()
    sns.boxplot(data=neut_sim, x=group_var, y="count", hue="pred_type")
    plt.title(group_var)
    plt.savefig("{}_neut_{}.svg".format(save_path, group_var))
    plt.show()
    plt.close()