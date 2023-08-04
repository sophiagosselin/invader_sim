use strict;
use warnings;
use File::Copy;
use Cwd qw(cwd);


my $directory = "ice_blast_runs";
my $threads = $ARGV[0];
system("rm -r unified_samples_experiment");

MAIN();

sub MAIN{
  mkdir("unified_samples_experiment");
  mkdir("unified_samples_experiment/1");
  open(my $allseqs, "+> unified_samples_experiment/1/all_sample_sequences.fasta");
  open(my $allinteins, "+> unified_samples_experiment/1/all_intein_sequences.fasta");

  opendir DIR, $directory ;
  while(my $dir = readdir DIR){
    # $dir will be the sample number
    if(-d "$directory\/$dir"){
      opendir DIR_SUB, "$directory\/$dir";
      while (my $dir_sub = readdir DIR_SUB){
        # $dir_sub will be the parameter directory
        if(-d "$directory\/$dir\/$dir_sub"){
          next;
        }
        elsif($dir_sub eq "invaded_seqs.fasta"){
          my %invaded_sample = READIN_FASTA("$directory\/$dir\/$dir_sub");
          foreach my $key (keys %invaded_sample){
            print $allseqs "$key\_sample_$dir\n";
            print $allseqs "$invaded_sample{$key}\n";
          }
        }
        elsif($dir_sub=~/intein_sample_.*/){
          my %intein_sequences = READIN_FASTA("$directory\/$dir\/$dir_sub");
          #my($file_handle)=($dir_sub=~/(.*?)\.fasta/);
          foreach my $key (keys %intein_sequences){
            print $allinteins "$key\_sample_$dir\n";
            print $allinteins "$intein_sequences{$key}\n";
          }
        }
        else{
          next;
        }
      }
    }
    else{
      next;
    }
  }

  close $allseqs;
  close $allinteins;

  my(%all_inteins)=READIN_FASTA("unified_samples_experiment/1/all_intein_sequences.fasta");
  my($random_intein_fasta)=RANDOM_INTEIN(\%all_inteins);
  my($intein_subset_fasta)=INTEIN_SUBSET($random_intein_fasta,\%all_inteins);

  my($int_subset_db)=MAKE_BLAST_DATABASE($intein_subset_fasta);
  my($extein_db)=MAKE_BLAST_DATABASE("unified_samples_experiment/1/all_sample_sequences.fasta");

  my @evals=("1e-3","1e-5","1e-8","1e-10","1e-13","1e-15","1e-18","1e-20");
  my @ids=(".5",".55",".6",".65",".7",".75",".8",".85",".9",".95");

  my $home_dir = cwd;
  foreach my $evalue (@evals){
    foreach my $identity (@ids){
      my($id_for_dir)=($identity=~/\.(.*)/);
      my $e_for_dir=$evalue;
      $e_for_dir=~s/\-/\_/g;
      my $directory = "$id_for_dir\_$e_for_dir";
      mkdir("unified_samples_experiment\/1\/$directory");
      #print "$dir\n";

      opendir RUN_DIR, "unified_samples_experiment/1";
      while(my $file = readdir RUN_DIR){
        if(-d "unified_samples_experiment/1/$file"){
          next;
        }
        else{
          #print "$file\n";
          copy("unified_samples_experiment/1/$file","$home_dir\/unified_samples_experiment\/1\/$directory\/$file");
        }
      }
      copy("iceblast.pl","unified_samples_experiment\/1\/$directory\/iceblast.pl");
      chdir "$home_dir\/unified_samples_experiment\/1\/$directory";
      system("export BLASTDB=$home_dir\/unified_samples_experiment\/1\/$directory");
      ICEBLAST("random_intein.fasta","all_sample_sequences.fasta","intein_subset.fasta",$identity,$evalue);
      chdir "$home_dir";
    }
  }
}

sub ICEBLAST{
  my $query_file = shift;
  my $search_db = shift;
  my $training_db = shift;
  my $id_param = shift;
  my $e_param = shift;
  system("perl iceblast.pl -in $query_file -psidb $training_db -outdb $search_db -t $threads -id $id_param -e $e_param -ds .25 -v 0");
}

sub INTEIN_SUBSET{
  my $int_fasta = shift;
  my %inteins = %{ my $hashref = shift};
  my %inteins_from_same_file;
  my $sample_of_int;

  my(%int_of_int)=READIN_FASTA($int_fasta);
  foreach my $asc (keys %int_of_int){
    ($sample_of_int)=($asc=~/.*?(intein_sample\_.*)/);
    #print "Random asc is $asc\t$sample_of_int\n";
  }

  foreach my $asc (keys %inteins){
    if($asc=~/$sample_of_int/){
      $inteins_from_same_file{$asc}=$inteins{$asc};
    }
    else{
      next;
    }
  }

  my @ascs_of_interest = (keys %inteins_from_same_file);
  my %random_subset;

  for(my $rng_counter=0; $rng_counter<10; $rng_counter++){
    my $rand_index = int rand @ascs_of_interest;
    my $rand_key = $ascs_of_interest[$rand_index];
    #print "Random index is $rand_index. Associated element is $rand_key\n";
    $random_subset{$rand_key}=$inteins_from_same_file{$rand_key};
    splice(@ascs_of_interest,$rand_index,1);
	}

  open(my $subset, "+> unified_samples_experiment\/1\/intein_subset.fasta");
  foreach my $print_asc (keys %random_subset){
    print $subset "$print_asc\n";
    print $subset "$random_subset{$print_asc}\n";
  }
  close $subset;

  return("unified_samples_experiment\/1\/intein_subset.fasta");
}

sub MAKE_BLAST_DATABASE{
  my $input_file = shift;
  system("makeblastdb -in $input_file -dbtype prot -parse_seqids -out $input_file");
  return("$input_file");
}

sub RANDOM_INTEIN{
  my %intein_seqs = %{my $hashref = shift};

  my @int_ascs = ( keys %intein_seqs );

  my $random_asc = $int_ascs[rand @int_ascs];

  open(my $rand_int, "+> unified_samples_experiment/1/random_intein.fasta");
  print $rand_int "$random_asc\n";
  print $rand_int "$intein_seqs{$random_asc}\n";
  close $rand_int;

  return("unified_samples_experiment/1/random_intein.fasta");
}

sub READIN_FASTA{
  #takes name of fasta formatted file as input
  #returns hash of sequences with ascessions as keys
  my $filehandle = shift;
  my %seq_data;
  my $asc_holder="";
  open(IN, "< $filehandle");
  while(<IN>){
    chomp;
    if($_=~/\>/){
      $asc_holder=$_;
      $seq_data{$asc_holder}="";
    }
    else{
      $seq_data{$asc_holder}.=$_;
    }
  }
  close IN;
  return(%seq_data);
}
