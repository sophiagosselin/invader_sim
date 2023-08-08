#!/usr/bin/env bash
#SBATCH --job-name=sim_seqs_subsample
#SBATCH --nodes=1
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=12gb
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin@uconn.edu
#SBATCH -o sim_%j.out
#SBATCH -e sim_%j.err

echo hostname
#NOTES: Need ice_blast.pl in home directory of simulation run
#as well as invader_sim.pl and a iSs.param file


#dependencies
module load perl/5.36.0
module load paml

#add libraries to path
export PERL5LIB=/home/FCAM/sgosselin/perl5/lib/perl5

#invader simulator. REMEMBER is.param file needs to be present
perl invader_sim.pl
