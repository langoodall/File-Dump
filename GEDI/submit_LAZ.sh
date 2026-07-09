#!/bin/bash
#SBATCH --job-name=LiDAR_Download	# Name of the job in the queue
#SBATCH --partition=msismall		# The queue I am putting it in
#SBATCH --time=48:00:00			# Maximum run time (HH:MM:SS)
#SBATCH --cpus-per-task=1		# Number of CPU cores requested
#SBATCH --mem=2G			# Memory requested
#SBATCH --output=output_%j.out		# Output file
#SBATCH --error=error_%j.err		# Error file

echo "Job started on $(date)"

# Move to the project directory
cd /projects/standard/freli001/shared/Louis_Projects/GEDI/submit

# Download all files listed in the laz_urls.txt
wget \
	-i laz_urls.txt \
	-P ../LAZ \
	-w 5 \
	-c \
	--tries 10 \
	--timeout 30

echo "Job ended on $(date)"
