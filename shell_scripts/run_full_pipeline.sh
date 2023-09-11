#!/usr/bin/env bash
#SBATCH --job-name=simulation_pipe
#SBATCH --nodes=1
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin@uconn.edu
#SBATCH -o sim_%j.out
#SBATCH -e sim_%j.err

#parameters for CPU and MEM based off requriements from full pipeline with experimental parameters in the manuscript

echo hostname
#NOTES: Need ice_blast.pl in home directory of simulation run
#as well as invader_sim.pl

#dependencies
module load perl/5.36.0
module load paml
module load blast/2.11.0
module load usearch

#add libraries to path
export PERL5LIB=/home/FCAM/sgosselin/perl5/lib/perl5

perl complete_pipeline.pl 4
