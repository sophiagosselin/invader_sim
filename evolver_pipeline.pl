#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;
use Getopt::Long;

#NEXT STEPS
#I might need to estimate parameters from a tree using MrBayes...
#Add a process for jumping the intein between exteins that isn't totally random...
#maybe use a normal distribution to determine the chances of jumping between a current branch and a nearby one? 
#BUT Make it so that every once in a while  there is a chance it will jump very far away 

#To make this work, will need the following:
#-Function to determine the distanc between two tips on a tree
#-Function to determine where the intein will jump to based on a probability distribution influenced by the patristic distances
#-Subroutine that, starting from 1 sequence progressively allows the intein to jump to new nodes on the tree 1 -> 2 -> 4.... 50 until the desired number is hit

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
  rename($extein_evolver_tree, "extein_tree.tre");
  $extein_evolver_tree = "extein_tree.tre";
  my $pdf_tree = WRITE_TREE($extein_evolver_tree);
}


sub EXTEIN_TREE{
  #make the extein tree
  open(my $evolver_pipe, '|-', 'paml-evolver') or die "Could not open pipe to paml-evolver: $!\n";
  print $evolver_pipe "1\n"; #uses random unrooted tree option
  print $evolver_pipe "$species_number\n"; #species number in tree
  print $evolver_pipe "$tree_number $random_seed\n"; #number of trees and random seed
  print $evolver_pipe "1\n"; #includes branch length from birth-death process
  print $evolver_pipe "$birth_rate $death_rate $sample_fraction $mutation_rate\n"; #see variable name
  print $evolver_pipe "0\n";
  close $evolver_pipe;
  #output will be evolver.out
  return("evolver.out");
}

sub WRITE_TREE{
  #takes newick tree as input
  #outputs tree as a pdf
  my $newick_file = shift;
  open(my $r_pipe, '|-', 'R --no-save') or die "Could not open pipe to R: $!\n";
  print $r_pipe "library(ape)\n";
  print $r_pipe "newick <- ape::read.tree(\"$newick_file\")\n";
  print $r_pipe "pdf(file=\"$newick_file.pdf\")\n";
  print $r_pipe "plot(newick)\n";
  print $r_pipe "dev.off()\n";
  print $r_pipe "quit(save = \"no\")\n";
  close $r_pipe;
  #returns pdf name
  return("$newick_file.pdf");
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
