#!/usr/bin/env bash
#SBATCH --job-name=ice_simulation
#SBATCH --nodes=1
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH -t 100:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin@uconn.edu
#SBATCH -o ice_sim_%j.out
#SBATCH -e ice_sim_%j.err

echo hostname
#NOTES: Need ice_blast.pl in home directory of simulation run
#as well as invader_sim.pl and a iSs.param file


#dependencies
module load blast
module load perl
module load R
module load paml

#add libraries to path
export PERL5LIB=/home/FCAM/sgosselin/perl5/lib/perl5

#invader simulator. REMEMBER is.param file needs to be present
perl invader_sim.pl

#get all simulation files from simulation
invaded_sequences=(~/invaded_sequences/*)
intein_sequences=(~/simulated_sequences/intein/*/*.fasta)
extein_sequences=(~/simulated_sequences/extein/*/*.fasta)

#create directory for ICE BLAST and prep variables
mkdir "ice_blast_runs"
counter=0
sample_directories=()
intein_files=()
extein_files=()

for file in $invaded_sequences
do
  #create subdirectory
  mkdir $counter

  #invaded_sequences/extein\_$file_descriptor1\_intein\_$file_descriptor2.fasta
  (extein,intein)=[[$file =~ invaded_sequences/extein\_(.*?)\_intein\_(.*?)\.fasta]]

  #move files to new subdirectory
  mv $file "ice_blast_runs/$counter/"
  mv "/simulated_sequences/intein/*/$intein.fasta" "ice_blast_runs/$counter/$intein.fasta.intein"
  mv "/simulated_sequences/extein/*/$extein.fasta" "ice_blast_runs/$counter/$extein.fasta.extein"
  cp "ice_blast.pl" "ice_blast_runs/$counter/"

  #push to array of subdirectory and inteins
  sample_directories+=("ice_blast_runs/$counter/")
  intein_files+=("$counter/$intein.fasta.intein")
  extein_files+=("$counter/$extein.fasta.extein")

  #increment up
  $counter++
done

#create a fasta file containing a random intein sequence from each sample

for intein_file in $intein_files
do
  #objects to hold file info
  declare -A paired_intein_sequneces
  intein_asc=()
  asc_holder=""
  (file_loc)=[[$intein_file =~ (.*)\/.*?\.fasta.intein]]

  #read fasta file to memory
  while read -r line;
    do
    if [[ $line=~^\> ]]
    then
      intein_asc+=($line)
      asc_holder=$line
    else
      paired_intein_sequneces[$asc_holder]+=$line
    fi

    #select random intein
    random_intein=${intein_asc[ $RANDOM % ${#intein_asc[@]} ]}

    #print random intein to file with asc
    echo $random_intein >> "$file_loc/random_intein.fasta"
    echo $paired_intein_sequneces[$random_intein] >> "$file_loc/random_intein.fasta"

    #now create subset of data for psidb
    #this value will need editing. Base it off the num of inteins used in dataset
    for rand in {0..50}
    do
      random_intein=${intein_asc[ $RANDOM % ${#intein_asc[@]} ]}
      #print random intein to file with asc
      echo $random_intein >> "$file_loc/intein_subset.fasta"
      echo $paired_intein_sequneces[$random_intein] >> "$file_loc/intein_subset.fasta"

      #the following ensures no repeats
      delete=($random_intein)
      asc_holder=( "${asc_holder[@]/$delete}")
    done
  done < $intein_file

  #finally, make blast db
  cd $file_loc
  makeblastdb -in $intein_file -input_type "prot" -parse_seqids -out "intein_sub.db"
  cd ..
  cd ..
done

#make blast databases from all extein sequences
for extein_file in $extein_files
do
  (file_loc)=[[$intein_file =~ (.*)\/.*?\.fasta.extein]]
  cd $file_loc
  makeblastdb -in $extein_file -input_type "prot" -parse_seqids -out "extein.db"
  cd ..
  cd ..
done

#create paired parameters to run ice_blast with
evals=("1e-2.5","1e-5","1e-7.5","1e-10","1e-12.5","1e-15","1e-17.5","1e-20")
paired_parameters=()
for id in $(seq .5 .55 .6 .65 .7 .75 .8 .85 .9 .95)
do
  for e in $evals
  do
    paired_parameters+="$id\&$e"
  done
done

#run ice_blast
for directory in $sample_directories
do
  directory_contents=(~/$directory/*)
  cd $directory
  for parameters in $paired_parameters
  do
    mkdir $parameters
    for file in $directory_contents
    do
      cp $file $paramaeters
    done

    cd $parameters
    (id_param)=[[$parmeters =~ (.*?)\&.* ]]
    (e_param)=[[$parmeters =~ .*?\&(.*) ]]

    perl iceblast.pl -in "random_intein.fasta" -psidb "intein_sub.db" -outdb "extein.db" -t 16 -id $id_param -e $e_param -ds .25

    cd ..
  done
  cd ..
  cd ..
done
