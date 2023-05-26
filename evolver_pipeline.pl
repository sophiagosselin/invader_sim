#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;
use Getopt::Long;
use Bio::TreeIO

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
my($species_number,$tree_number,$birth_rate,$death_rate,$sample_fraction,$mutation_rate,%inputs);
my $input_primer = "Parameters should be in the following format:\n
sp - # of species.
tn - # of trees.
b - birth rate.
d - death rate.
sf - sample fraction.
m - mutation rate.\n
Example: sp 50 b 1 d 1 sf .1 m .01\n\n";

#Code Start

MAIN();

sub MAIN{
  #starts by creating the extein tree.
  my $user_check = "N";
  my $tree_type = "extein";
  my $input_primer_append = "Please input parameters for the $tree_type phylogeny now!\n".$input_primer;
  %inputs = parse_and_check_inputs($input_primer_append);
  #%inputs applies to the specific run. Therefore critical to keep seperate
  #begin parameter testing for extein phylogeny
  print "Beginning extein simulation parameters test.\n Testing will conclude once the user is happy with the tree and set parameters.\n"
  tree_simulation_parameter_testing($user_check,$tree_type);
  #save parameters settled on
  my %extein_inputs =  %inputs;
  #now evolves the intein tree
  #gets new parameters first
  $tree_type = "intein";
  $input_primer_append = "Please input parameters for the $tree_type phylogeny now!\n".$input_primer;
  %inputs = parse_and_check_inputs($input_primer_append);
  #begin parameter testing for intein phylogeny
  print "Beginning intein simulation parameters test.\n Testing will conclude once the user is happy with the tree and set parameters.\n"
  tree_simulation_parameter_testing($user_check,$tree_type);
  #save parameters settled on
  my %extein_inputs =  %inputs;
  #once the user is satisfied with all parameters, ask for number of experiments to run
  print "Parameters have been set; ready to run experiment.\n How many samples would you like to create?\n\n"
  my $number_of_samples = <STDIN>;
  chomp $number_of_samples;
  #next up loop tree creating and intein insertion process
}

sub parse_and_check_inputs{
  #takes a string as the message for the user and primes them to input options
  #then splits apart that response, and assigns values to global variables
  my $message = shift;
  print($message);
  my $new_options = <STDIN>;
  chomp $new_options;
  my @new_options_list = split /\s+/, $new_options;
  GetOptions ('sp=s' => \$species_number, 'b=s' =>\$birth_rate,
              'd=s' => \$death_rate, 'sf=s' => \$sample_fraction,
              'm=s' => \$mutation_rate, @new_options_list);
  #Dicitonary for checking input parameters
  my %inputs_for_tree = ("species_number" => $species_number,
                "tree_number" => $tree_number,
                "birth_rate" => $birth_rate,
                "death_rate" => $death_rate,
                "sample_fraction" => $sample_fraction,
                "muation_rate" => $mutation_rate);
    #check input variables
  foreach my $input (keys %inputs_for_tree){
    if(!defined $inputs_for_tree{$input}){
      print "$input is not defined. Please define it now:";
      my $new_parameter = <STDIN>;
      chomp $new_parameter;
      $inputs_for_tree{$input}=$new_parameter;
    }
    else{}
  }
  return(%inputs_for_tree);
}

sub tree_simulation_parameter_testing{
  #takes a y/n toggle, tree extsension and prossibly tree object as inputs
  #loops through tree building process until the user is satisfied based on the
  #global parameters
  #returns nothing, but the final chosen parameters will be in the global hash %inputs
  my $toggle = shift;
  my $tree_extension = shift;
  my $hashref = shift || "";
  if(uc($toggle) eq "N"){
    my $evolver_tree = tree_simulation();
    rename($evolver_tree, "$tree_extension"."_tree.tre");
    $evolver_tree = "$tree_extension"."_tree.tre";
    my $pdf_tree = write_tree($evolver_tree);
    print "Does this simulated $tree_extension tree look good to you [Y/N]?\n";
    $toggle = <STDIN>;
    chomp $toggle;
    if(uc($toggle) eq "N"){
      print "Would you like to change your tree building parameters [Y/N]?\n";
      my $new_param_toggle = <STDIN>;
      chomp $new_param_toggle;
      if(uc($new_param_toggle) eq "Y"){
        %inputs = parse_and_check_inputs($input_primer);
      }
      $hashref= tree_simulation_parameter_testing($toggle,$tree_extension);
    }
    elsif(uc($toggle) eq "Y"){
      return();
    }
    else{
      print "$toggle is not Y or N. Please print Y or N to continue...\n";
      $toggle = <STDIN>;
      chomp $toggle;
      tree_simulation_parameter_testing($toggle,$tree_extension,$tree_to_return);
    }
  }
  elsif(uc($toggle) eq "Y"){
    return();
  }
  else{
    print "$toggle is not Y or N. Please print Y or N to continue...\n";
    $toggle = <STDIN>;
    chomp $toggle;
    $hashref = tree_simulation_parameter_testing($toggle,$tree_extension,$tree_to_return);
  }
  return();
}

sub tree_simulation{
  #simulated a tree based on parameters provided
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

sub write_tree{
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

sub show_pdf {
  #takes pdf file as input. Then displays PDF to screen
  my $pdf_file = shift;
  my $pdf_txt;
  open (my $fh, $pdf_file);
  # set the file handle to binary mode
  binmode $fh;
  # read it all into a string;
  while (<$fh>){
    $pdf_txt .= $_;
  }
  close ($fh);

  #takes string and displays PDF
  my $method = "Content-disposition:inline; filename='$pdf_file'"; # default method
  my $size = length($pdf);
  print "Content-Type: application/pdf\n";
  print "Content-Length: $size\n";
  print "$method\n";
  print "Content-Transfer-Encoding: binary\n\n"; # blank line to separate headers
  print $pdf;
}


sub readin_newick{
  #readin newick format tree file
  #return a tree object and an array of tips
  my $tree_file = shift;
  my $tree_input = new Bio::TreeIO(-file => $tree_file,
                            -format => "newick");
  my $tree_object = $tree_input->next_tree;
  my @tree_taxa = $tree_object->get_leaf_nodes;
  return($tree_object,@tree_taxa);
}

sub random_array_element{
  #readin array
  #return random element from said array
  my @array = @_;
  my $random_entry = $array[rand @array];
  return($random_entry);
}

sub get_patristic_distance{
  #takes two tree nodes and the associated tree object as inputs
  #returns the distance between said nodes
  my $node1 = shift;
  my $node2 = shift;
  my $tree_object = shift;
  my $distance = $tree_object->distance(-nodes => [$node1, $node2]);
  return($distance);
}

sub find_nearest_neighbors{
  #takes an array of tree tip nodes, a node of interest, and a tree object as inputs
  #returns a hash of tips paired to their distance from the node of interest
  my $tree_object = shift;
  my $tip_of_interest = shift;
  my @array_of_tips = @_;
  my %distances;
  foreach my $tip (@array_of_tips){
    if($tip eq $tip_of_interest){
      next;
    }
    else{
      my $node_distance = get_patristic_distance($tip,$tip_of_interest);
      $distances{$tip}=$node_distance;
    }
  }
  return(%distances);
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
