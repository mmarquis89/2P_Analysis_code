#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 01:20:00                   		# Runtime in minutes
#SBATCH -p priority                			# Partition (queue) to submit to
#SBATCH --mem=2G               				# memory needed

expDate=$1
sid=$2

echo make_fictrac_vids $expDate $sid

module load matlab/2017b
matlab -nodesktop -r "make_fictrac_vids('$expDate', $sid)"
