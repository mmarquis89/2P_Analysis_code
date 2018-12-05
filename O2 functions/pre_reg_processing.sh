#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 6:00:00                   		# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=200G
               			# memory needed
parentDir=$1
sid=$2
expDate=$3
imgSaveDir=$4

echo pre_reg_processing $expDate $sid $parentDir $imgSaveDir



module load matlab/2017a
matlab -nodesktop -r "preReg_routine_MM('$parentDir', '$sid', '$expDate', 'OutputDir', '$imgSaveDir')"