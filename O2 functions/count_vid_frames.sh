#!/bin/bash

#SBATCH -c 1                   			    # Number of cores requested
#SBATCH -t 1:00:00                   			# Runtime in minutes
#SBATCH --mem=2G
#SBATCH -p short                			# Partition (queue) to submit to

echo count_vid_frames
vidDir=$1
sid=$2

module load matlab/2017a
myVar=$(matlab -nodesktop -r "count_vid_frames('$vidDir', $sid)")
