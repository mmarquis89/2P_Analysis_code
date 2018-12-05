#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 00:10:00                   		# Runtime in minutes
#SBATCH -p short            			    # Partition (queue) to submit to
#SBATCH --mem=2G             			    # memory needed (memory PER CORE)

expDate=$1
sid=$2

module load matlab/2017a
matlab -nodesktop -r "O2_preprocessing_test('$expDate', $sid)"

