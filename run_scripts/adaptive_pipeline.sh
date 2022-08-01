#!/bin/bash

while getopts t:r:u:d:o: flag; do
  case "$flag" in
    t) TREE_FILE=${OPTARG};;
    r) SINGLE_REGIME_FILE=${OPTARG};;
    u) TWO_REGIME_FILE=${OPTARG};;
    d) DATA_PATH=${OPTARG};;
    o) OUTPUT_FILE=${OPTARG};;
  esac
done

Rscript run_scripts/adaptive_evogenex.R $TREE_FILE $SINGLE_REGIME_FILE $TWO_REGIME_FILE $DATA_PATH $OUTPUT_FILE
