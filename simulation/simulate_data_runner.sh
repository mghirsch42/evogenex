#!/bin/bash

tree_file="tree_files/resolved/sc-bwes-cons-resolved-10.tree"
outpath="test.csv"
regime_file="regime_files/resolved/agg_resolved.csv"
theta_root=150
theta_ratios=( 1.06 1.04 1.02 )
sigmasqs=( .15 1.5 15 )
gammas=( .05 .5 1 )
alphas=( .125 1 8 )
ngene=2
nrep=8
rep_drop=.1

for theta_ratio in "${theta_ratios[@]}"; do
    theta_small=$(echo "(2*$theta_root)/(1+$theta_ratio)" | bc -l)
    theta_large=$(echo "(2*$theta_root*$theta_ratio) / (1+$theta_ratio)" | bc -l)
    for sigmasq in "${sigmasqs[@]}"; do
        for gamma in "${gammas[@]}"; do
            for alpha in "${alphas[@]}"; do
                outfile="${outpath}tr_${theta_ratio}_sq_${sigmasq}_g_${gamma}_a_${alpha}.csv"
                echo $theta_ratio " " $sigmasq " " $gamma " " $alpha
                Rscript simulate_data.R \
                    $tree_file \
                    $outfile \
                    $regime_file \
                    $theta_root \
                    $theta_small \
                    $theta_large \
                    $sigmasq \
                    $gamma \
                    $alpha \
                    $ngene \
                    $nrep \
                    $rep_drop \
                    "nonaggressive" \
                    "aggressive"
                exit
            done
        done
    done
done