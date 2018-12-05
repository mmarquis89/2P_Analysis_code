#!/bin/sh

#SBATCH -c 1                   			    # Number of cores requested
#SBATCH -t 2:00:00                  		# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=50G               			# memory needed (memory PER CORE)
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  

parentDir=$1
sid=$2
imgSaveDir=$3

echo make_fluorescence_vid $parentDir $sid $imgSaveDir

module load matlab/2017a
matlab -nodesktop -r "create_average_fluorescence_vid('$parentDir', '$sid', 'OutputDir', '$imgSaveDir')"