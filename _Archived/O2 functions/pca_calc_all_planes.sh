#!/bin/bash

#SBATCH -c 1                   			    # Number of cores requested
#SBATCH -t 3-00:00:00                		# Runtime in minutes
#SBATCH -p medium                			# Partition (queue) to submit to
#SBATCH --mem=200G               			# memory needed (memory PER CORE)

parentDir=$1
sessionDataFile=$2

echo pca_calc $parentDir $sessionDataFile

module load matlab/2017a

matlab -nodesktop -r "pca_calc_all_planes('$parentDir', '$sessionDataFile')"
#matlab -nodesktop -r "pca_calc('$parentDir', '$sessionDataFile')"