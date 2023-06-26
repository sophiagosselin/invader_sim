#!/usr/bin/env bash
#SBATCH --job-name=sim_test
#SBATCH --nodes=1
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH -t 100:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin@uconn.edu
#SBATCH -o sim_%j.out
#SBATCH -e sim_%j.err

echo hostname

#dependencies
module load perl/5.36.0
module load R
module load paml

#add libraries to path
#this leads to bioperl libraries
export PERL5LIB=/home/FCAM/sgosselin/perl5/lib/perl5

#invader simulator. REMEMBER is.param file needs to be present
perl invader_sim.pl
