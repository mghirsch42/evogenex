#!/bin/bash

while getopts t:r:d:o: flag; do
  case "$flag" in
    t) TREE_FILE=${OPTARG};;
    r) REGIME_FILE=${OPTARG};;
    d) DATA_PATH=${OPTARG};;
    o) OUTPUT_FILE=${OPTARG};;
  esac
done

Rscript run_evogenex.R $TREE_FILE $REGIME_FILE $DATA_PATH $OUTPUT_FILE