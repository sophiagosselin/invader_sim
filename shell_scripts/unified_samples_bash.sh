#!/usr/bin/env bash
#SBATCH --job-name=ice_all_samples
#SBATCH --nodes=1
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20gb
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin@uconn.edu
#SBATCH -o all_samples_sim_%j.out
#SBATCH -e all_samples_sim_%j.err

echo hostname

#dependencies
module load blast/2.11.0
module load perl
module load usearch

perl unified_samples_experiment.pl 16
