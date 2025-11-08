#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=7-00:00:00
#SBATCH --mem=48G

####### Set environment variables ###############
#
start_time=$(date +%s)
echo "Job started on: $(date)"
####### Run your script #########################
eval "$(conda shell.bash hook)"
conda activate guideseq

GUIDEseq_DIR="/work/vetmed_data/jj/projects/basilHubbard/guideseq/guideseq/guideseq"
python ${GUIDEseq_DIR}/guideseq.py main -m manifest.yaml
#seqkit seq --reverse --complement guideseq_S0_R2_001.fastq -o guideseq_S0_R2_001_rc.fastq

end_time=$(date +%s)
echo "Job ended on: $(date)"

# Calculate the difference in times to get the duration
duration=$((end_time - start_time))

hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))
echo "Job duration: $hours hours, $minutes minutes, $seconds seconds"
