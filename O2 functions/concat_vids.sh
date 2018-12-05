#!/bin/bash

#SBATCH -c 1                   			    # Number of cores requested
#SBATCH -t 1:00:00                			# Runtime in minutes
#SBATCH -p priority                			# Partition (queue) to submit to
#SBATCH --mem=16G               			# memory needed (memory PER CORE)

vidDir=$1
sid=$2

echo concat_vids $vidDir $sid


fileStr="*sid_${sid}_tid*.avi"
outputFileName="sid_${sid}_AllTrials"

module load matlab/2017a
matlab -nodesktop -r "concat_vids('$vidDir', '$fileStr', 'OutputFile', '$outputFileName')"