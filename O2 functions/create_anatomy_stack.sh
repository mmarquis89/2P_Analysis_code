#!/bin/sh

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 2:00:00                   	    # Runtime in minutes
#SBATCH -p short               			    # Partition (queue) to submit to
#SBATCH --mem=32G               			# memory needed (memory PER CORE)


echo create_anatomy_stack

parentDir=$1
imgSaveDir=$2

echo "$parentDir $imgSaveDir" >> test.txt

module load matlab/2017a
matlab -nodesktop -r "create_anatomy_stack('$parentDir', 'OutputDir', '$imgSaveDir', 'FileString', '*Stack_*.tif')"