#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 06:00:00                   		# Runtime in minutes
#SBATCH -p priority              			# Partition (queue) to submit to
#SBATCH --mem=2G             			    # memory needed (memory PER CORE

vidDataDir=$1
vidSaveDir=$2
sid=$3

echo make_behavior_vids $vidDataDir $vidSaveDir $sid

module load matlab/2017a
matlab -nodesktop -r "make_behavior_vids('$vidDataDir', '$vidSaveDir', $sid)"

# Create text file with list of trials/tids
# module load matlab/2017a
# matlab -nodesktop -r "create_trial_list('$vidDataDir', $sid)"
# sleep 10

#------------ Start a job to create each video -------------------------------------------------
# file="$vidDataDir/sid_${sid}_dirList.txt"
# while IFS=',' read -r dirName tid
# do
	# sbatch make_vid.sh $dirName $vidSaveDir $sid $tid 
# done <"$file"
# sleep 5

# ------------ Concatenate vids afterwards -----------------------------------------------------
# jid1=$(sbatch --dependency=singleton --job-name=makeVid concat_vids.sh $vidSaveDir $sid)
# jid1="${jid1//[!0-9]/}"

# ---- Count the number of frames in each of the individual trial videos -----------------------
# jid2=$(sbatch --dependency=singleton --job-name=makeVid count_vid_frames.sh $vidSaveDir $sid)
# jid2="${jid2//[!0-9]/}"
# sleep 2

# ------------- Calculate optic flow for all trials --------------------------------------------
# file="$vidDataDir/sid_${sid}_dirList.txt"
# echo $vidSaveDir
# while IFS=',' read -r dirName tid
# do
	# echo $vidSaveDir $sid $tid
	# sbatch optic_flow_calc.sh $vidSaveDir $sid $tid
# done <"$file"
# sleep 10

# sbatch --dependency=singleton --job-name=makeVid optic_flow_calc.sh $vidSaveDir $sid $tid
# jid3=$(sbatch --dependency=singleton --job-name=makeVid optic_flow_calc.sh $vidSaveDir $sid)
# jid3="${jid3//[!0-9]/}"
# sleep 2

# ------------- Normalize optic flow data across all trials ------------------------------------
# jid2=$(sbatch normalize_optic_flow.sh $vidSaveDir $sid)
# jid2="${jid2//[!0-9]/}"

# ------------- Create optic flow vids for anvil annotation ------------------------------------
# file="$vidDataDir/sid_${sid}_dirList.txt"
# while IFS=',' read -r dirName tid
# do
	# sbatch --dependency=afterok:$jid2 create_optic_flow_vid.sh $vidSaveDir $sid $tid
# done <"$file"
# sleep 10

# ------------ Concatenate optic flow vids -----------------------------------------------------
# sbatch concat_optic_flow_vids.sh $vidSaveDir $sid 

# sbatch --dependency=singleton --job-name=makeVid concat_optic_flow_vids.sh $vidSaveDir $sid 




