import pandas as pd
import os
import argparse

# This will save a file with the following columns:
# 'ensemble_id', 'gene_name', 'ou1_conv', 'ou1_theta', 'ou1_alpha',
#       'ou1_sigma_sq', 'ou1_gamma', 'ou1_loglik', 'ou2_conv', 'ou2_theta',
#       'ou2_theta_base', 'theta_diff', 'ou2_alpha', 'ou2_sigma_sq',
#       'ou2_gamma', 'ou2_loglik', 'brown_conv', 'brown_theta',
#       'brown_sigma_sq', 'brown_gamma', 'brown_loglik', 'ou2_vs_bm_pvalue',
#       'ou2_vs_ou1_pvalue'

def main(result_path, save_path, model, inverse):
    # Set column names and key words based on if it's a constrained or adaptive experiment
    if model == "c":
        model_col = "constrained_vs_neutral"
        cols = ["ensemble_id", "gene_name", "pvalue", "qvalue", model_col]
        if inverse:
            model_goal = "neutral"
        else:
            model_goal = "constrained"
    if model == "a":
        model_col = "adaptive"
        cols = ["ensemble_id", "gene_name"]
        if inverse:
            model_goal = "not-adaptive"
        else:
            model_goal = "adaptive"

    results_df = pd.DataFrame(columns=cols)

    # Read all the data
    # Any columns we didn't define above will be added
    for f in [f for f in os.listdir(result_path) if ".csv" in f]:
        temp_df = pd.read_csv(result_path+f)
        results_df = pd.concat([results_df, temp_df], ignore_index=True)
        # results_df = results_df.append(temp_df, ignore_index=True)
    # The results will have the ensemble id and gene name as a single column - split into two
    # The simulation will not have this, just set each key to the gene name
    if results_df["gene"].astype(str).str.contains("_").all():
        results_df[["ensemble_id", "gene_name"]] = results_df["gene"].str.split("_", expand=True)
    else:
        results_df["ensemble_id"] = results_df["gene"]
        results_df["gene_name"] = results_df["gene"]
    results_df = results_df.drop("gene", axis=1)
    
    # If we have an adaptive experiment, it will have two rows, one with the ou2 base theta and one with the ou2 adaptive theta
    # So we need to combine those into one
    if model == "a":
        new_cols = results_df.columns.to_list()
        theta_idx = new_cols.index("ou2_theta")
        new_cols.insert(theta_idx+1, "ou2_theta_base")
        new_cols.insert(theta_idx+2, "theta_diff")
        df2 = pd.DataFrame(columns=new_cols)
        for g in results_df["gene_name"].unique():
            sub = results_df[results_df["gene_name"] == g]  # Two entries - adaptive, then base
            df2.loc[len(df2)] = sub.iloc[0]
            # df2 = df2.append(sub.iloc[0]) # Add adaptive entry
            df2.loc[df2["gene_name"] == g, "ou2_theta_base"] = sub["ou2_theta"].iloc[1] # Add base value
    
        # Calculate the difference between the theta values
        df2["theta_diff"] = df2["ou2_theta"] - df2["ou2_theta_base"]

        results_df = df2
    
    # Get only the adaptive results
    pos_res = results_df[results_df[model_col] == model_goal]
    pos_res = pos_res.drop(model_col, axis=1)

    print(str(len(pos_res)) + " positive results out of " + str(len(results_df)))

    if save_path:
        # Save all results and the positive results
        pos_res.to_csv(save_path + "gene_info.csv", index=False)
        results_df.to_csv(save_path + "all_results.csv", index=False)
    else:
        print(pos_res)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("result_path", type=str, action="store")
    parser.add_argument("model", type=str, action="store")
    parser.add_argument("-s", "--save_path", type=str, action="store", default=None)
    parser.add_argument("-i", "--inverse", action="store_true")
    args = parser.parse_args()
    if args.model not in ["c", "a"]:
        print("Model must be either 'c' for constrained or 'a' for adaptive.")
        exit()
    main(args.result_path, args.save_path, args.model, args.inverse)
