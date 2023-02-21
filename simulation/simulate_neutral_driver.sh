#!/bin/bash

tree_file=$1
outpath=$2

theta_root=150
#sigmasqs=( .15 1.5 15 )
sigmasqs=( .15 15 )
#rs=( 10 100 1000 )
rs=( 10 1000 )
ngene=100
nrep=5
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
