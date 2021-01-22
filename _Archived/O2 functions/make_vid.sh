#!/bin/bash

#SBATCH -c 1                   			    # Number of cores requested
#SBATCH -t 2:00:00                 			# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=4G               			    # memory needed (memory PER CORE)
#SBATCH --job-name=makeVid

vidDataDir=$1
vidSaveDir=$2
sid=$3
tid=$4

echo make_vid $vidDataDir $vidSaveDir $sid $tid

module load matlab/2017a
matlab -nodesktop -r "make_vid('$vidDataDir', $sid, $tid, 'OutputDir', '$vidSaveDir')"