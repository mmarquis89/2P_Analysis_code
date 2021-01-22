#!/bin/sh

#SBATCH -t 5                  				# Runtime in minutes
#SBATCH -p priority              			# Partition (queue) to submit to
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  

expDate=$1
sid=$2

echo preprocessing $expDate $sid

#expDate="2018_04_14_exp_1"
#sid=0

imgDataDir="/n/scratch2/mjm60/${expDate}/ImagingData"
vidDataDir="/n/scratch2/mjm60/${expDate}/BehaviorVideo"
imgSaveDir="/n/scratch2/mjm60/${expDate}/sid_${sid}/ImagingData"
vidSaveDir="/n/scratch2/mjm60/${expDate}/sid_${sid}/BehaviorVideo"

# Make average fluorescence vids
sbatch make_fluorescence_vid.sh "$imgDataDir" "$sid" "$imgSaveDir"

# Make anatomy stack
sbatch create_anatomy_stack.sh "$imgDataDir" "$imgSaveDir"

#-----------------------------------------------------------------
# Process Imaging Data
#-----------------------------------------------------------------

jid1=$(sbatch pre_reg_processing.sh $imgDataDir $sid $expDate $imgSaveDir)
jid1="${jid1//[!0-9]/}"

sbatch --dependency=afterok:$jid1 NoRMCorre_registration.sh $imgSaveDir $sid $expDate

#-----------------------------------------------------------------
# Create Behavior Vids
#-----------------------------------------------------------------

# Create behavior vids
sbatch make_behavior_vids.sh $vidDataDir $vidSaveDir $sid






