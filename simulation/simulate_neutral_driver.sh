#!/bin/bash

tree_file=$1
outpath=$2

theta_root=150
sigmasqs=( .15 1.5 15 )
rs=( 100 150 200 )
ngene=100
nrep=8
rep_drop=.1

for sigmasq in "${sigmasqs[@]}"; do
    for r  in "${rs[@]}"; do
        outfile="${outpath}sq_${sigmasq}_r_${r}.csv"
        echo $sigmasq " " $r
        Rscript simulation/simulate_neutral.R \
            $tree_file \
            $outfile \
            $theta_root \
            $sigmasq \
            $r \
            $ngene \
            $nrep \
            $rep_drop
    done
done
