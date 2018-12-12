#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 00:10:00                   		# Runtime in minutes
#SBATCH -p priority              			# Partition (queue) to submit to
#SBATCH --mem=8G             			    # memory needed 

expDate=$1
sid=$2

# echo analysis_processing $expDate $sid

# imgDataDir="/n/scratch2/mjm60/${expDate}/ImagingData"
# vidDataDir="/n/scratch2/mjm60/${expDate}/BehaviorVideo"
# imgSaveDir="/n/scratch2/mjm60/${expDate}/sid_${sid}/ImagingData"
# vidSaveDir="/n/scratch2/mjm60/${expDate}/sid_${sid}/BehaviorVideo"

# sessionDataFile="rigid_sid_${sid}_sessionFile.mat"
# echo $sessionDataFile

module load matlab/2017a
matlab -nodesktop -r "analysis_processing('$expDate', $sid)"



# ------------- Create and save an analysis metadata file ---------------
# jid1=$(sbatch initial_analysis_processing.sh $imgSaveDir $sessionDataFile)
# jid1="${jid1//[!0-9]/}"

# ------------------ Calculate overall behavior state dF/F ---------------------------
# sbatch  --dependency=afterok:$jid1 behavioral_state_dff_calc.sh $imgSaveDir $sessionDataFile

# ------------------ Calculate and save PCA data ---------------------------
# sbatch --dependency=afterok:$jid1 pca_calc.sh $imgSaveDir $sessionDataFile