#!/bin/bash

#SBATCH -c 6                    			# Number of cores requested
#SBATCH -t 8:00:00                   		# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=230G               			# memory needed


parentDir=$1
sid=$2

echo normcorre_registration $parentDir $sid

fileName="sid_${sid}_sessionFile.mat"

module load matlab/2017a
matlab -nodesktop -r "normcorre_registration('$parentDir', '$fileName')"