#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 00:05:00                			# Runtime in minutes
#SBATCH -p priority            			    # Partition (queue) to submit to
#SBATCH --mem=2G             			    # memory needed (memory PER CORE

expDate=$1
sid=$2

echo calc_event_dff $expDate $sid


imgSaveDir="/n/scratch2/mjm60/${expDate}/sid_${sid}/ImagingData"
fileStr="EventData*.mat"
sessionDataFile="rigid_sid_${sid}_sessionFile.mat"
echo $sessionDataFile

module load matlab/2017a
matlab -nodesktop -r "calc_event_dff_avg_startScript('$imgSaveDir', '$fileStr', '$sessionDataFile')"