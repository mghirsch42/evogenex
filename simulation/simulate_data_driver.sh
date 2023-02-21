#!/bin/bash

#tree_file="tree_files/resolved/sc-bwes-cons-resolved-10.tree"
#outpath="test.csv"
#regime_file="regime_files/resolved/agg_resolved.csv"

tree_file=$1
regime_file=$2
outpath=$3

theta_root=150
theta_ratios=( 1.1 1.08 1.06 1.04 1.02 )
#theta_ratios=( 1.1 1.06 1.02 )
#sigmasqs=( .15 1.5 15 )
sigmasqs=( 37.5 )
rs=( 10 100 1000 )
#rs=( 10 1000 )
alphas=( .125 1 8 16 )
#alphas=( 16 )
ngene=100
nrep=8
#nrep=5
rep_drop=.1

for theta_ratio in "${theta_ratios[@]}"; do
    theta_small=$(echo "(2*$theta_root)/(1+$theta_ratio)" | bc -l)
    theta_large=$(echo "(2*$theta_root*$theta_ratio) / (1+$theta_ratio)" | bc -l)
    for sigmasq in "${sigmasqs[@]}"; do
        for r  in "${rs[@]}"; do
            for alpha in "${alphas[@]}"; do
                outfile="${outpath}tr_${theta_ratio}_sq_${sigmasq}_r_${r}_a_${alpha}.csv"
                echo $theta_ratio " " $sigmasq " " $r " " $alpha
                Rscript simulation/simulate_data.R \
                    $tree_file \
                    $outfile \
                    $regime_file \
                    $theta_root \
                    $theta_small \
                    $theta_large \
                    $sigmasq \
                    $r \
                    $alpha \
                    $ngene \
                    $nrep \
                    $rep_drop \
                    "nonaggressive" \
                    "aggressive"
            done
        done
    done
done
