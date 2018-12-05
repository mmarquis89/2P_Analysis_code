#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 02:00:00                   		# Runtime in minutes
#SBATCH -p priority                			# Partition (queue) to submit to
#SBATCH --mem=200G               			# memory needed

expDate=$1
sid=$2

echo extract_ROI_data $expDate $sid

sessionDataFile="rigid_sid_${sid}_sessionFile.mat"
ROIfile="ROI_metadata.mat"
imgSaveDir="/n/scratch2/mjm60/${expDate}/sid_${sid}/ImagingData"

module load matlab/2017a
matlab -nodesktop -r "extract_ROI_data('$imgSaveDir', '$sessionDataFile', '$ROIfile')"
