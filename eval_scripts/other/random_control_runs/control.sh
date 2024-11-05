
# {x..y} includes y
for i in {1..50}; do
    python eval_scripts/other/random_control_runs/random_regime_file.py -n $i
    mkdir results/control/random${i}/
    mkdir results/control/gene_lists/random${i}
    ./run_scripts/adaptive_runner.sh \
        -t tree_files/sc-bwes-cons-resolved-10.tree \
        -r regime_files/control/single.csv \
        -u regime_files/control/random${i}.csv \
        -d data/tpm/ \
        -o control/random${i}/ \
        -s
done