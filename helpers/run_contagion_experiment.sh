#!/bin/bash

# Script to run the non-linear contagion
# experiment for each network in 'datasets',
# for each sampler in 'samplers'

# Define an array of strings
datasets=("lyon" "high" "email-Eu" "email-Enron")
samplers=("NuDHy")

# Iterate through the array
for samp in "${samplers[@]}"; do
	for db in "${datasets[@]}"; do
		echo "Running GENEPY and ECI"
		echo "---- `date`"
		python run_non_linear_contagion.py $db $samp
	done
done