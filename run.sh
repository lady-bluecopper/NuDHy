#!/usr/bin/env bash

echo ''
echo ''
echo '  ______                                ______           '
echo ' |  ____|                              |  ____|          '
echo ' | |____    _____    __  __    _____   | |____    _____  '
echo ' |____  |  |____ |  |  \/  |  |  _  |  |____  |  |____ | '
echo '  ____| |   / __ |  | |\/| |  | |_| |   ____| |   / __ | '
echo ' |______|  |_____|  |_|  |_|  |_____|  |______|  |_____| '
echo '|""""""""||"""""""||""""""""||"""""""||""""""""||"""""""|'
echo -e '\n\n'

# Loading configurations for experiments
echo '>> Loading config file config.cfg'

source config.cfg

unset datasets
declare -A datasets
datasets[$ecoli]=$ecoli_defaults
datasets[$iaf]=$iaf_defaults
datasets[$ijo]=$ijo_defaults
datasets[$enron]=$enron_defaults
datasets[$ord]=$ord_defaults
datasets[$hs]=$hs_defaults
datasets[$dblp9]=$dblp9_defaults
datasets[$cit]=$cit_defaults
datasets[$math]=$math_defaults

unset flags
declare -A flags
flags[$ecoli]=$ecoli_flags
flags[$iaf]=$iaf_flags
flags[$ijo]=$ijo_flags
flags[$enron]=$enron_flags
flags[$ord]=$ord_flags
flags[$hs]=$hs_flags
flags[$dblp9]=$dblp9_flags
flags[$cit]=$cit_flags
flags[$math]=$math_flags

unset threshs
declare -A threshs
threshs[$ecoli]=$ecoli_threshs
threshs[$iaf]=$iaf_threshs
threshs[$ijo]=$ijo_threshs
threshs[$enron]=$enron_threshs
threshs[$ord]=$ord_threshs
threshs[$hs]=$hs_threshs
threshs[$dblp9]=$dblp9_threshs
threshs[$cit]=$cit_threshs
threshs[$math]=$math_threshs

echo -e '\n\n>> Creating directories ...'
mkdir -p $resultsDir

for dataset in ${!datasets[@]}
do
	datasetPath=${datasetsDir}/${dataset}
	default=${datasets[${dataset}]}
	flag=${flags[${dataset}]}
	thresh=${threshs[${dataset}]}
	defaults=(`echo $default|tr "," "\n"`)
	experiments=(`echo $flag|tr "," "\n"`)

	echo ">> Processing dataset ${dataset} with default values (${defaults[@]})"
	echo ">> Experiment flags ${experiments[@]}"

	if [[ ${experiments[0]} -eq "1" ]]; then
		echo '-----------------------'
		echo '      Convergence      '
		echo '-----------------------'

		OUTPUT="$resultsDir/convergence/"
		mkdir -p $OUTPUT

		echo "Running command ..."
		echo "$JVM $CONV_jar datasetsDir=$datasetsDir datasetName=$dataset resultsDir=$resultsDir seed=$seed minFreqs=$thresh k=$k size=$size"
		echo "---- `date`"
		$JVM $CONV_jar datasetsDir=$datasetsDir datasetName=$dataset resultsDir=$resultsDir seed=$seed minFreqs=$thresh k=$k size=$size
	fi


	if [[ ${experiments[1]} -eq "1" ]]; then
		echo '-----------------------------------'
		echo '        Structural Entropy         '
		echo '-----------------------------------'

		OUTPUT="$resultsDir/entropy/"
		mkdir -p $OUTPUT

		echo "Running command ..."
		echo "$JVM $ENT_jar datasetsDir=$datasetsDir datasetName=$dataset resultsDir=$resultsDir numThreads=$numThreads seed=$seed numSwaps=${defaults[0]} numSamples=${defaults[1]} samplerType=$samplerType"
		echo "---- `date`"
		$JVM $ENT_jar datasetsDir=$datasetsDir datasetName=$dataset resultsDir=$resultsDir numThreads=$numThreads seed=$seed numSwaps=${defaults[0]} numSamples=${defaults[1]} samplerType=$samplerType
	fi

        if [[ ${experiments[2]} -eq "1" ]]; then
                echo '-----------------------------------'
                echo '            Scalability            '
                echo '-----------------------------------'

                OUTPUT="$resultsDir/scalability/"
                mkdir -p $OUTPUT

                echo "Running command ..."
                echo "$JVM $SCALA_jar datasetsDir=$datasetsDir datasetName=$dataset resultsDir=$resultsDir numSwaps=${defaults[0]} numSamples=${defaults[1]} samplerType=$samplerType"
                echo "---- `date`"
                $JVM $SCALA_jar datasetsDir=$datasetsDir datasetName=$dataset resultsDir=$resultsDir numSwaps=${defaults[0]} numSamples=${defaults[1]} samplerType=$samplerType
        fi

        if [[ ${experiments[3]} -eq "1" ]]; then
                echo '-----------------------------------'
                echo '             Sampling              '
                echo '-----------------------------------'

                OUTPUT="$resultsDir/sampling/"
                mkdir -p $OUTPUT

                echo "Running command ..."
                echo "$JVM $SAMP_jar datasetsDir=$datasetsDir datasetName=$dataset resultsDir=$resultsDir maxNumSwapsFactor=$maxNumSwapsFactor numSamples=${defaults[1]} samplerType=$samplerType headTailDisjoint=$headTailDisjoint"
                echo "---- `date`"
                $JVM $SAMP_jar datasetsDir=$datasetsDir datasetName=$dataset resultsDir=$resultsDir maxNumSwapsFactor=$maxNumSwapsFactor numSamples=${defaults[1]} samplerType=$samplerType headTailDisjoint=$headTailDisjoint
        fi
done
echo 'Terminated.'
