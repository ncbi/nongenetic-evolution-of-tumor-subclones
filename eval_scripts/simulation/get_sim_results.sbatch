#!/bin/bash

#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=36gb
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=output/%x_%j.out
#SBATCH --error=output/%x_%j.err

module load Python3/3.8.1

source ../env/bin/activate

RESULT_PATH=$1

for f in $RESULT_PATH*; do
    # Create folders for each parameter set result if they don't exist
    if [ ! -d $f ]; then
        exp_folder="${f%.csv}"
        mkdir $exp_folder
        mv $f $exp_folder
    else
        exp_folder=$f
    fi
    # If we don't already have the adaptive results collated, run the script for that
    if [ ! -f ${exp_folder}gene_info.csv ]; then
        echo $exp_folder
        python3 eval_scripts/get_results/get_gene_results.py $exp_folder/ a -s $exp_folder/
    fi
done

# Note: must load R after running Python script
module load R/4.1.2/3.14

# Add stats columns to the results files
for exp_folder in $RESULT_PATH*/; do
    if [ -f ${exp_folder}gene_info.csv ]; then 
        echo $exp_folder
        Rscript eval_scripts/get_results/add_stats_cols.R ${exp_folder}gene_info.csv
        Rscript eval_scripts/get_results/add_stats_cols.R ${exp_folder}all_results.csv
    fi
done

deactivate

