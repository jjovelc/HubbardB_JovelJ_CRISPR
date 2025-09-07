#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=7-00:00:00
#SBATCH --mem=16G

####### Handle user-defined prefix #########################

# Ensure the user provides a prefix
if [ -z "$1" ]; then
  echo "Usage: $0 PREFIX"
  exit 1
fi

PREFIX=$1
LOG_DIR="${PREFIX}_log-err_files"
mkdir -p "$LOG_DIR"

# Redirect output and error logs
#SBATCH --output=%x.out
#SBATCH --error=%x.err

exec > "$LOG_DIR/${PREFIX}.log"
exec 2> "$LOG_DIR/${PREFIX}.err"

####### Set environment variables ###############

start_time=$(date +%s)
echo "Job started on: $(date)"

####### Run your script #########################
eval "$(conda shell.bash hook)"
conda activate cas9_nickase

bash run_cutadapt.sh "$PREFIX"

end_time=$(date +%s)
echo "Job ended on: $(date)"

# Calculate the difference in times to get the duration
duration=$((end_time - start_time))

hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))
echo "Job duration: $hours hours, $minutes minutes, $seconds seconds"

