#!/bin/bash

#SBATCH -c 1                   			    # Number of cores requested
#SBATCH -t 5:00:00                			# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=16G               			# memory needed (memory PER CORE)

vidDir=$1
sid=$2

flowVidDir="$vidDir/opticFlowVids"
echo concat_optic_flow_vids $vidDir $sid $flowVidDir

fileStr="*sid_${sid}_tid*_With_Optic_Flow.avi"
outputFileName="sid_${sid}_AllTrials_With_Optic_Flow"

module load matlab/2017a
matlab -nodesktop -r "concat_vids('$flowVidDir', '$fileStr', 'OutputFile', '$outputFileName')"