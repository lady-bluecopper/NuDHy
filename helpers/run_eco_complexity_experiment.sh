#!/bin/bash

# Script to compute GENEPY and ECI
# for each trade network in 'datasets',
# for each sampler in 'samplers'

# Define an array of strings
datasets=("hs1995_thresh" "hs2009_thresh" "hs2019_thresh" "hs2020_thresh")
samplers=("NuDHy")

# Iterate through the array
for samp in "${samplers[@]}"; do
	for db in "${datasets[@]}"; do
		echo "Running GENEPY and ECI"
		echo "---- `date`"
		python run_genepy.py $db $samp
	done
done
