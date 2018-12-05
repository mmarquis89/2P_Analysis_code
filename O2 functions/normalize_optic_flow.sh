#!/bin/bash

#SBATCH -c 1                   			    # Number of cores requested
#SBATCH -t 00:30:00                   		# Runtime in minutes
#SBATCH --mem=1G
#SBATCH -p priority                			# Partition (queue) to submit to


vidSaveDir=$1
sid=$2

echo normalize_optic_flow vidSaveDir $vidSaveDir sid $sid

module load matlab/2017a
myVar=$(matlab -nodesktop -r "normalize_optic_flow('$vidSaveDir', $sid)")