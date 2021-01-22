#!/bin/bash

#SBATCH -c 1                   			    # Number of cores requested
#SBATCH -t 11:00:00                   		# Runtime in minutes
#SBATCH --mem=2G
#SBATCH -p short                			# Partition (queue) to submit to

echo optic_flow_calc
vidDir=$1
sid=$2
tid=$3

roiDataFile="Behavior_Vid_ROI_Data.mat"

echo "optic flow - $vidDir $roiDataFile $sid $tid"

module load matlab/2017a
myVar=$(matlab -nodesktop -r "single_trial_optic_flow_calc('$vidDir', $sid, $tid, '$roiDataFile')")