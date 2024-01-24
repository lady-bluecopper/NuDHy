#!/bin/bash

# Scripts to compute each metric for each dataset in 'datasets'
# and for each sampler in 'samplers'

# Array of datasets of interest (without the extension)
datasets=("eco01100" "metabolic_iaf1260b" "metabolic_ijo1366" "uspto19762016" "dblp_v9" "citation_software" "enron" "qna_math")
# Array of samplers of interest
samplers=("NuDHy" "Base" "BaseD" "UnpretentiousNullModel")

for samp in "${samplers[@]}"; do
	for db in "${datasets[@]}"; do
		echo "Running Reciprocity"
		echo "---- `date`"
		$1 run_reciprocity.py $db $samp

		echo "Running Shell Indices Computation"
		echo "---- `date`"
		$1 run_coreness.py $db $samp

		echo "Running Multi-Order Laplacian"
		echo "---- `date`"
		$1 run_laplacian.py $db $samp

		echo "Running Centrality"
		echo "---- `date`"
		$1 run_centrality.py $db $samp

		echo "Running Hyper-coreness Computation"
		echo "---- `date`"
		$1 compute_hypercoreness.py $db $samp
	done
done
