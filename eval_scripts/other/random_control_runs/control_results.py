import pandas as pd

##########
## Compile the control results into a single file
##########

result_path = "results/control/gene_lists/random{}/gene_info.csv"
regime_path = "regime_files/control/random{}.csv"

df = pd.DataFrame()

for i in range(1, 51):
    print(i)
    regime = pd.read_csv(regime_path.format(i))
    chosen_nodes = regime[regime["regime"]=="1_chosen"]["node"].tolist()
    results = pd.read_csv(result_path.format(i))
    l.append([i] + chosen_nodes + [len(results), len(results[results["logFC"]>0]), len(results[results["logFC"]<0])])

df = pd.DataFrame(data=l, columns=["run", "node1", "node2", "node3", "total", "up", "down"])
df.to_csv("results/control/results.csv".format(n), index=False)
