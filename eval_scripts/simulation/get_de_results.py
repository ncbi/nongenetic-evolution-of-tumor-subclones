import pandas as pd
import os
import argparse

def main(result_path, save_path):
    results_df = pd.DataFrame()
    sim = result_path.split("/")[-2]

    # Read all the data
    # Any columns we didn't define above will be added
    for f in [f for f in os.listdir(result_path) if ".csv" in f]:
        temp_df = pd.read_csv(result_path+f)
        temp_df["sim"] = f[:-4]
        results_df = results_df.append(temp_df, ignore_index=True)
    results_df = results_df.drop("ensemble_id", axis=1)
    print(results_df)
    results_df["de"] = (results_df["fdr_p"] < 0.05).astype(int)
    # count = results_df.groupby(sim_vars).count()
    count = results_df.groupby("sim").sum().reset_index()
    count = count[["sim", "de"]]
    print(count)

    sim_params = count["sim"].str.split("_")
    print(sim_params)
    if sim == "neut_sim":
        count["sq"] = sim_params.str[1]
        count["r"] = sim_params.str[3]
        count = count[["sq", "r", "de"]]
    if sim == "const_sim":
        count["sq"] = sim_params.str[3]
        count["r"] = sim_params.str[5]
        count["a"] = sim_params.str[7]
        count = count[["sq", "r", "a", "de"]]
    if sim == "adpt_sim":
        count["theta_ratio"] = sim_params.str[1]
        count["sq"] = sim_params.str[3]
        count["r"] = sim_params.str[5]
        count["a"] = sim_params.str[7]
        count = count[["sq", "r", "a", "theta_ratio", "de"]]
    count.to_csv(save_path+sim+"_de.csv", index=False)
    print(count)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--result_path", type=str, action="store", default="results/simulated/differential_expression/neut_sim/")
    parser.add_argument("-s", "--save_path", type=str, action="store", default="results/simulated/summaries/")
    args = parser.parse_args()

    main(args.result_path, args.save_path)
