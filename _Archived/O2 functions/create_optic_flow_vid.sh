#!/bin/bash

#SBATCH -c 1                   			    # Number of cores requested
#SBATCH -t 11:00:00                			# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=128G             			    # memory needed

vidSaveDir=$1
sid=$2
tid=$3

outputDir="$vidSaveDir/opticFlowVids"

echo create_optic_flow_vid $vidSaveDir $sid $tid $outputDir

frameCountFile="sid_${sid}_frameCounts.txt"

module load matlab/2017a
matlab -nodesktop -r "create_single_trial_optic_flow_vid('$vidSaveDir', $sid, $tid, '$frameCountFile', 'OutputDir', '$outputDir')"
