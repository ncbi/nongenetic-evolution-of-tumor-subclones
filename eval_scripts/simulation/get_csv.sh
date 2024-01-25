#!/bin/bash

results_dir="results/simulated/"
result_file="results/simulated/all_counts.csv"

echo sim_type,tree,regime,n_rep,theta_ratio,sq,nb_r,alpha,count > $result_file

for sim_dir in $results_dir*/; do
    sim="${sim_dir%/}"
    sim="${sim##*/}"
    regime_tree_repl=(${sim//_/ })
    echo $sim | grep "neut"
    if echo "$sim" | grep -q "neut"; then
        sim_type="neut"
        regime=${regime_tree_repl[2]}
        tree=${regime_tree_repl[3]}
        repl=${regime_tree_repl[4]}
    else
        sim_type="adpt"
        regime=${regime_tree_repl[1]}
        tree=${regime_tree_repl[2]}
        repl=${regime_tree_repl[3]}
    fi
    tree="${tree#t}"
    repl="${repl#r}"
    echo $regime $tree $repl
    for exp_dir in $sim_dir*/; do
        exp="${exp_dir%/}"
        exp="${exp##*/}"
        tr_sq_r_a=(${exp//_/ })
        if [ "${sim_type}" = "neut" ]; then
            tr=""
            sq=${tr_sq_r_a[1]}
            nbr=${tr_sq_r_a[3]}
            alpha=""
        else
            tr=${tr_sq_r_a[1]}
            sq=${tr_sq_r_a[3]}
            nbr=${tr_sq_r_a[5]}
            alpha=${tr_sq_r_a[7]}
        fi
        n=$(wc -l < ${exp_dir}gene_info.csv)
        n=$((n-1))
        # n=$(awk -F, 'BEGIN { sum = 0 } $10 < 0.05 { ++sum } END { print sum }' ${exp_dir}de.csv)
        # echo "${exp_dir}de.csv"
        awk -F, '{ print $10 }' "${exp_dir}de.csv"
        echo $n
        #exit
        echo "${sim_type},${tree},${regime},${repl},${tr},${sq},${nbr},${alpha},${n}" >> $result_file
    done
done
