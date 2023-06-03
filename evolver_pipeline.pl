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
  #get user inputs and test until they are satisfied
  my %extein_inputs = test_parameters("extein");
  my %intein_inputs = test_parameters("intein");

  #once the user is satisfied with all parameters, ask for number of experiments to run
  print "Parameters have been set; ready to run experiment.\n How many samples would you like to create?\n\n"
  my $number_of_samples = <STDIN>;
  chomp $number_of_samples;

  #create the number of trees for inteins and exteins based on the number of experiments requested
  my @extein_phylogenies = simulate_n_trees("extein",$number_of_samples,%extein_inputs);
  my @intein_phylogenies = simulate_n_trees("intein",$number_of_samples,%intein_inputs);

  #pair each extein tree with a random intein tree. Each pairing constitutes a sample
  my %paired_intein_extein_trees;
  foreach my $extein_tree (@extein_phylogenies){
    my $random_intein_tree = splice(@intein_phylogenies,@intein_phylogenies[rand($intein_phylogenies)],1);
    $paired_intein_extein_trees{$extein_tree}=$random_intein_tree;
  }

  #foreach paired sample invade the exteins with the inteins
  invade_exteins(%paired_intein_extein_trees);
}

sub test_parameters{
  #takes tree type as input, returns user determined inputs for that tree type.
  my $tree_type = shift;
  my $user_check = "N";
  my $input_primer_append = "Please input parameters for the $tree_type phylogeny now!\n"."$input_primer";
  %inputs = parse_and_check_inputs($input_primer_append);
  #begin parameter testing for phylogeny of interest
  print "Beginning extein simulation parameters test.\n Testing will conclude once the user is happy with the tree and set parameters.\n"
  parameter_testing_loop($user_check,$tree_type);
  #save parameters settled on
  my %downstream_inputs =  %inputs;
  return(%downstream_inputs)
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

sub parameter_testing_loop{
  #takes a y/n toggle, tree extsension and prossibly tree object as inputs
  #loops through tree building process until the user is satisfied based on the
  #global parameters
  #returns nothing, but the final chosen parameters will be in the global hash %inputs
  my $toggle = shift;
  my $tree_extension = shift;
  if(uc($toggle) eq "N"){
    my $evolver_tree = simulate_1_tree();
    my $pdf_tree = write_tree($evolver_tree);
    print "Does this simulated $tree_extension tree look good to you [Y/N]?\n";
    $toggle = <STDIN>;
    chomp $toggle;
    if(uc($toggle) eq "N"){
      %inputs = parse_and_check_inputs($input_primer);
      parameter_testing_loop($toggle,$tree_extension);
    }
    elsif(uc($toggle) eq "Y"){
      return();
    }
  }
  elsif(uc($toggle) eq "Y"){
    return();
  }
  return();
}

sub simulate_n_trees{
  #takes hash of evolver inputs, a file handle, and the number of simulations as input
  #returns array of tree files
  my $file_handle = shift;
  my $simulation_num = shift;
  my %inputs_for_tree = shift;
  my @tree_files;
  #make directory to store trees in
  mkdir($file_handle);

  #change global inputs to ones needed
  %inputs = %inputs_for_tree;
  my @loop = (1..$simulation_num);
  foreach my $num (@loop){
    my $phylogeney = simulate_1_tree();
    my $path_to_file = "$file_handle\$file_handle\_$num\.tree";
    copy($phylogeney,$path_to_file);
    push(@tree_files,$path_to_file);
    unlink $phylogeney;
  }
  return(@tree_files);
}

sub simulate_1_tree{
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

sub invade_exteins{
  #takes a set of paired intein and extein trees as input
  #invades the extein trees with the inteins using a monte carlo chain
  #returns a nested hash for which the key is the extein tree, and has
  #data for which tips are invaded, and with which intein as nested data
  my %tree_samples = shift;
  foreach my $extein_tree (keys %tree_samples){
    my($extein_tree_object,@extein_tips)=readin_newick($extein_tree);
    my($intein_tree_object,@intein_tips)=readin_newick($tree_samples{$extein_tree});

    #get pairwise distance between all nodes in intein and extein tree
    my %extein_pairwise_distances= pairwise_patristic_distance($extein_tree_object,@extein_tips);
    my %intein_pairwise_distances= pairwise_patristic_distance($intein_tree_object,@intein_tips);

    #select a random intein and random extein to invade w/ said intein
    my $intein_tip = @intein_tips[rand($intein_tips)];
    my $extein_tip = @extein_tips[rand($extein_tips)];

    #take the random intein/extein pair as the starting point and begin the MC chain
    #test if there is still an intein that hasn't invaded a sequence.
    my %paired_sequences;
    until(!@intein_tips){
      $paired_sequences{$extein_tip}=$intein_tip;
      my $closest_intein = get_smallest_distance($intein_tip,%intein_pairwise_distances);
      my %distance_from_current_extein = get_distances_from_tip($extein_tip,%extein_pairwise_distances);
      ### Keep going here!







      
    }
  }
}

sub get_distances_from_tip{
  #takes a nested hash of distances and key1
  #the hash is formatted as key1->key2->value
  #returns a hash of key2's with the associated distance to key1
  my $key1 = shift;
  my %distance_hash = shift;
  foreach my $key (keys %distance_hash){
    next if($key ne $key1);
    my %return_hash;
    foreach my $key2 (keys %{$distance_hash{$key1}}){
      next if($key2 eq $key);
      $return_hash{$key2}=$distance_hash{$key1}->{$key2};
    }
  }
  return(%return_hash);
}

sub get_smallest_distance{
  #takes a nested hash of distances and key1
  #the hash is formatted as key1->key2->value
  #returns the key2 with the smallest distance to key1
  my $key1 = shift;
  my %distance_hash = shift;
  foreach my $key (keys %distance_hash){
    next if($key ne $key1);
    my $smallest = "toggle";
    my $key_holder;
    foreach my $key2 (keys %{$distance_hash{$key1}}){
      if(!$smallest="toggle"){
        if($smallest < $distance_hash{$key1}->{$key2}){
          next;
        }
      }
      $smallest = $distance_hash{$key1}->{$key2};
      $key_holder = $key2;
    }
  }
  return($key_holder);
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

sub get_patristic_distance{
  #takes two tree nodes and the associated tree object as inputs
  #returns the distance between said nodes
  my $node1 = shift;
  my $node2 = shift;
  my $tree_object = shift;
  my $distance = $tree_object->distance(-nodes => [$node1, $node2]);
  return($distance);
}

sub pairwise_patristic_distance{
  #takes an array of tree tip nodes, and a tree object as inputs
  #returns a nested hash of tips paired to their distance from the node of interest
  my $tree_object = shift;
  my @array_of_tips = @_;
  my %distances;
  foreach my $tip (@array_of_tips){
    foreach my $tip_to_compare (@array_of_tips){
      my $node_distance = get_patristic_distance($tip,$tip_to_compare);
      $distances{$tip}={$tip_to_compare}=$node_distance;
    }
  }
  return(%distances);
}

sub monte_carlo_chain{
  #takes a window width
  #and an array with associated distances from the current position
  #returns the new location
  my $window_width = shift;
  my @distance_of_tips = shift;
  my $current_position = 0;

  #uniform random deviant
  my $uniform_random = rand();

  #creates a integer value based on the window size
  #then gets new position proposal
  my $new_position = $current_position + int($uniform_random*$window_width);

  my $acceptor = accept_or_deny_move(@distance_of_tips[$new_position]);

  if($acceptor == 1){
    return($new_position);
  }
  else{
    return(0);
  }
}

sub accept_or_deny_move{
  #given an edge length for a move as input, find the probability of making that proposed position the new position
  #using an exponential distribution and a uniform random deviant
  #returns a boolean 1/0

  my $lenght_to_new_position = shift;
  #rate paramaeter for the exponential dist.
  my $lambda = 0.5;

  #prob of acceptance based on dist.
  my $probability_of_acceptance = ($lambda*(exp(-$lambda*$lenght_to_new_position)));

  my $uniform_random_deviant = rand();
  #compare uniform random to prob of acceptance. If lower, new position is accpeted!
  if($uniform_random_deviant <= $probability_of_acceptance){
    return(1);
  }
  else{
    return(0);
  }
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
