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
#NOTES: Need ice_blast.pl in home directory of simulation run
#as well as invader_sim.pl and a iSs.param file


#dependencies
module load blast
module load perl/5.36.0
module load R
module load paml

#add libraries to path
export PERL5LIB=/home/FCAM/sgosselin/perl5/lib/perl5

#invader simulator. REMEMBER is.param file needs to be present
perl invader_sim.pl

#get all simulation files from simulation
invaded_sequences=(invaded_sequences/*)

#create directory for ICE BLAST and prep variables
mkdir "ice_blast_runs"
counter=0
sample_directories=()
intein_files=()
extein_files=()

for file in ${invaded_sequences[@]}
do
  #create subdirectory
  mkdir "ice_blast_runs/$counter"

  #invaded_sequences/extein\_$file_descriptor1\_intein\_$file_descriptor2.fasta
  [[ $file =~ invaded_sequences\/extein_sample_(.*)_intein_sample_(.*)\.fasta ]]
  extein=${BASH_REMATCH[1]}
  intein=${BASH_REMATCH[2]}

  #move files to new subdirectory
  mv $file "ice_blast_runs/$counter/"
  mv  simulated_sequences/extein/$extein/*.fasta "ice_blast_runs/$counter/extein_sample_$extein.fasta"
  mv  simulated_sequences/intein/$intein/*.fasta "ice_blast_runs/$counter/intein_sample_$intein.fasta"
  cp iceblast.pl "ice_blast_runs/$counter/"

  #push to array of subdirectory and inteins
  sample_directories+=("ice_blast_runs/$counter")
  intein_files+=("ice_blast_runs/$counter/intein_sample_$intein.fasta")
  extein_files+=("ice_blast_runs/$counter/extein_sample_$extein.fasta")

  #increment up
  ((counter+=1))
done

#create a fasta file containing a random intein sequence from each sample

for intein_file in ${intein_files[@]}
do
  #objects to hold file info
  declare -A paired_intein_sequences
  intein_asc=()
  asc_holder=""
  [[ $intein_file =~ (.*)\/.*\.fasta ]]
  file_loc=${BASH_REMATCH[1]}
  sed -i -e 's/\///g' $intein_file
  sed -i -e 's/\-//g' $intein_file
  #read fasta file to memory
  count_2=0
  while read -r line;
  do
    line=$(echo $line)
    if [[ $line =~ ^\> ]]
    then
      ${asc_holder}="${line}_${count_2}"
      paired_intein_sequences[${asc_holder}]=""  # Use quotes around the index
      intein_asc+=(${asc_holder})
      count_2=$((count_2 + 1))
      echo "Valid ASC ID: ${asc_holder}"
    else
      subseq=${paired_intein_sequences[${asc_holder}]}
      seq="${subseq}${line}"
      paired_intein_sequences[${asc_holder}]="${seq}"  # Use quotes around the index
      #echo "Readin report: ${seq}"
    fi
  done < "$intein_file"

  #select random intein
  random_index=$[$RANDOM % ${#intein_asc[@]}]
  random_intein=${intein_asc[$random_index]}

  #print random intein to file with asc
  echo "${random_intein}" >> "$file_loc/random_intein.fasta"
  echo "Random asc: ${random_intein}"
  echo "Random seq:"
  rnd_seq=${paired_intein_sequneces[${random_intein}]}
  echo "${rnd_seq}" >> "$file_loc/random_intein.fasta"
  echo "${rnd_seq}"

  #now create subset of data for psidb
  #this value will need editing. Base it off the num of inteins used in dataset
  size_of_arr=${#intein_asc[@]}
  rand_size=$((size_of_arr / 10))
  for rand in {0..$rand_size}
    do
      random_index=${intein_asc[ $RANDOM % ${#intein_asc[@]} ]}
      random_intein=$intein_asc[$random_index]
      #print random intein to file with asc
      echo "${random_intein}" >> "${file_loc}/intein_subset.fasta"
      echo "Acession: ${random_intein}"
      seq_from_array=${paired_intein_sequences[$random_intein]}
      echo "${seq_from_array}" >> "${file_loc}/intein_subset.fasta"
      echo "Sequence: ${seq_from_array}"
      #the following ensures no repeats
      delete=($random_intein)
      asc_holder=( "${asc_holder[@]/$delete}")
    done

  #finally, make blast db
  makeblastdb -in "${file_loc}/intein_subset.fasta" -dbtype "prot" -parse_seqids -out "${file_loc}/intein_subset.fasta"

done

#make blast databases from all extein sequences
#first need to make each asc unique.
for extein_file in ${extein_files[@]}
do
  sed -i -e 's/\///g' $extein_file
  sed -i -e 's/\-//g' $extein_file
  count_3=0
  while read -r line;
  do
    if [[ $line =~ ^\> ]]
    then
      echo "${line}_${count_3}" >> "${extein_file}.temp"
      count_3=$((count_3 + 1))

    else
      echo "${line}" >> "${extein_file}.temp"
    fi
  done < "$extein_file"

  #rename and remove old fasta file
  rm $extein_file
  mv "${extein_file}.temp" $extein_file

  #now for db
  [[ $extein_file =~ (.*)\/(.*\.fasta) ]]
  file_loc=${BASH_REMATCH[1]}
  out_name=${BASH_REMATCH[2]}

  makeblastdb -in $extein_file -dbtype "prot" -parse_seqids -out "${file_loc}/all_exteins.fasta"
  mv $extein_file "${file_loc}/all_exteins.fasta"
done

#create paired parameters to run ice_blast with
evals=("1e-3" "1e-5" "1e-8" "1e-10" "1e-13" "1e-15" "1e-18" "1e-20")
paired_parameters=()
ids=(".5" ".55" ".6" ".65" ".7" ".75" ".8" ".85" ".9" ".95")
for id in ${ids[@]}
do
  for e in ${evals[@]}
  do
    paired_parameters+=("${id}&${e}")
  done
done

#run ice_blast
for directory in ${sample_directories[@]}
  do
    for parameters in ${paired_parameters[@]}
      do
        [[ $parameters =~ (.*)\&(.*) ]]
        id_param=${BASH_REMATCH[1]}
        e_param=${BASH_REMATCH[2]}
        paramdir=$parameters
        paramdir1=$(echo "${paramdir//\-/}")
        paramdir2=$(echo "${paramdir1//\./}")
        paramdir3=$(echo "${paramdir2//\&/}")
        cd $directory
        mkdir $paramdir3
        for file in *;
          do
            if [[ -f $file ]]
              then
                cp $file $paramdir3
              fi
          done

        cd $paramdir3
        perl iceblast.pl -in random_intein.fasta -psidb intein_subset.fasta -outdb all_exteins.fasta -t 8 -id $id_param -e $e_param -ds .25 -v 2
        wait
        cd -

        cd ..
        cd ..
      done
  done
