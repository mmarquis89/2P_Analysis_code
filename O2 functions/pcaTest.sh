#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 07:05:00                   		# Runtime in minutes
#SBATCH -p priority            			    # Partition (queue) to submit to
#SBATCH --mem=100G             			    # memory needed (memory PER CORE)


module load matlab/2017a
matlab -nodesktop -r "pcaTest()"

