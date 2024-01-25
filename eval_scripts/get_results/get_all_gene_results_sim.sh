#/bin/bash

##########
### Loop over all result files for a specific simulation and run
### get_gene_results.py and add_stats_cols.R
##########

result_path="results/simulated/adpt_sim_adpt_run/"

for f in $result_path*/; do
    echo $f
    python3 eval_scripts/get_results/get_gene_results.py $f a -s $f
    Rscript eval_scripts/get_results/add_stats_cols.R ${f}gene_info.csv
    Rscript eval_scripts/get_results/add_stats_cols.R ${f}/all_results.csv
done
