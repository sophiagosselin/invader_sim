#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;
use Getopt::Long;

#USAGE: species_num tree_num birth_rate death_rate sampling_fraction mutation_rate

#GLOBALS
my $random_seed = int rand(1000000000);
my($species_number,$tree_number,$birth_rate,$death_rate,$sample_fraction,$mutation_rate);
GetOptions ('sp=s' => \$species_number, 'tn=s' =>\$tree_number, 'b=s' =>\$birth_rate, 'd=s' => \$death_rate, 'sf=s' => \$sample_fraction, 'm=s' => \$mutation_rate);

#Dicitonary for checking input parameters
my %inputs = ("species_number" => $species_number,
              "tree_number" => $tree_number,
              "birth_rate" => $birth_rate,
              "death_rate" => $death_rate,
              "sample_fraction" => $sample_fraction,
              "muation_rate" => $mutation_rate);

#Code Start
SETUP(\%inputs);
MAIN();

sub SETUP{
  #check input variables
  my $hash_ref = shift;
  my %inputs_to_check = %{$hash_ref};
  foreach my $input (keys %inputs_to_check){
    if(!defined $inputs_to_check{$input}){
      die "$input is not defined. Please check your inputs!\n";
    }
  }
}

sub MAIN{
  my $extein_evolver_tree = EXTEIN_TREE();
  die;
}


sub EXTEIN_TREE{
  #make the extein tree
  system("paml-evolver");
  system("1"); #uses random unrooted tree option
  system("$species_number"); #species number in tree
  system("$tree_number\ $random_seed"); #number of trees and random seed
  system("1"); #includes branch length from birth-death process
  system("$birth_rate\ $death_rate\ $sample_fraction\ $mutation_rate"); #see variable name
  #output will be evolver.out

}


sub old_code{
my $counter=0;
open(IN, "< $ARGV[0]");
my $tree_string;
while(<IN>){
  chomp;
  $tree_string.=$_;
}
close IN;
my @trees=split(/\;/,$tree_string);
foreach my $tree (@trees){
  print "$tree\n";
  mkdir("$counter");
  my $random_seed = int(rand(100000));
  if(0 == $random_seed % 2){
    $random_seed++;
  }
  else{}
  open(EVO, "+> $counter/MCaa.dat");
  print EVO
"0\n
$random_seed\n
\n
5000 150 1\n
\n
-1\n
\n
$tree\;\n
\n
.5 4\n
2 /usr/lib/paml/data/dat/wag.dat\n
\n
0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05\n
0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05\n
\n
A R N D C Q E G H I\n
L K M F P S T W Y V\n
\n
\/\/ end of file";
  close EVO;
  chdir("$counter");
  system("paml-evolver 7");
  chdir("..");
  $counter++;
}
}
