#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH -c 1                 #optional: number of cpus, default is 1
#SBATCH --mem=50GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --mail-user=jkhay@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=dedup              #optional: job name
#SBATCH --output=dedup_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=dedup_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID



samtools sort -O sam /projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam -o C1_SE_uniqAlign_sorted.sam

/usr/bin/time -v ./Hays_deduper.py -u STL96.txt -f C1_SE_uniqAlign_sorted.sam -o dedup_out.sam

exit