#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;
use Getopt::Long;
use Bio::TreeIO;
use Bio::AlignIO;

#USAGE: perl evolver_pipeline.pl

#GLOBALS
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
  print "Parameters have been set; ready to run experiment.\n How many samples would you like to create?\n\n";
  my $number_of_samples = <STDIN>;
  chomp $number_of_samples;

  #create the number of trees for inteins and exteins based on the number of experiments requested
  my @extein_phylogenies = simulate_n_trees("extein",$number_of_samples,\%extein_inputs);
  my @intein_phylogenies = simulate_n_trees("intein",$number_of_samples,\%intein_inputs);

  #pair each extein tree with a random intein tree. Each pairing constitutes a sample
  my %paired_intein_extein_trees;
  foreach my $extein_tree (@extein_phylogenies){
    my $random_intein_tree = splice(@intein_phylogenies,@intein_phylogenies[rand($#intein_phylogenies)],1);
    $paired_intein_extein_trees{$extein_tree}=$random_intein_tree;
  }

  #ask the user for max jump size when invading from one extein tip to another
  print "Preparing intein invasion simulation.
  What window size would you like to constrain intein invasions to?
  Importantly, in an extein tree with n tips, and x inteins, this window should be at most n-x wide.\n
  For example, in a tree with 10 tips, a window of 5 would allow the intein to invade any of the closest 5 tips to its current position.\n";
  my $window_size = <STDIN>;
  chomp $window_size;

  #simulates intein invasion
  #returns a nested hash in the following format:
  #$extein tree -> {"tips"/"intein_tree"} -> {$extein_tips/$intein_file_name} -> $invading_intein
  my %paired_invasion_data = invade_exteins($window_size,\%paired_intein_extein_trees);

  #now simulate sequence evolution based on the phylogenies
  #returns a hash of the sequence file associated with the tree file
  my ($extein_size,%paired_extein_tree_and_sequence_files) = simulate_sequences_from_tree("extein",$extein_inputs{"species_number"},@extein_phylogenies);
  my ($intein_size,%paired_intein_tree_and_sequence_files) = simulate_sequences_from_tree("intein",$intein_inputs{"species_number"},@intein_phylogenies);

  #now invade extein sequences with intein sequences
  insert_sequences($extein_size,\%paired_invasion_data,\%paired_extein_tree_and_sequence_files,\%paired_intein_tree_and_sequence_files);

  print "Simulation finished! Find the fasta files of invaded extein sequences with their un-invaded counterparts in the 'invaded_sequences' directory.\n";
}

sub test_parameters{
  #takes tree type as input, returns user determined inputs for that tree type.
  my $tree_type = shift;
  my $user_check = "N";
  my $input_primer_append = "Please input parameters for the $tree_type phylogeny now!
  Remember that the extein tree should have at least 2x the amount of taxa as the intein tree.
  (This is critical for the method in which this software simulates intein invasion)\n\n"."$input_primer";
  %inputs = parse_and_check_inputs($input_primer_append);
  #begin parameter testing for phylogeny of interest
  print "Beginning extein simulation parameters test.\n Testing will conclude once the user is happy with the tree and set parameters.\n";
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
  my %inputs_for_tree = %{my $hashref = shift};
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
  #first generates a random seed and ensures it is an odd value
  my $random_tree_seed = int rand(1000000000);
  if(0 == $random_tree_seed % 2){
    $random_tree_seed++;
  }
  open(my $evolver_pipe, '|-', 'paml-evolver') or die "Could not open pipe to paml-evolver: $!\n";
  print $evolver_pipe "1\n"; #uses random unrooted tree option
  print $evolver_pipe "$species_number\n"; #species number in tree
  print $evolver_pipe "$tree_number $random_tree_seed\n"; #number of trees and random seed
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
  my $pdf_text;
  open (my $fh, $pdf_file);
  # set the file handle to binary mode
  binmode $fh;
  # read it all into a string;
  while (<$fh>){
    $pdf_text .= $_;
  }
  close ($fh);

  #takes string and displays PDF
  my $method = "Content-disposition:inline; filename='$pdf_file'"; # default method
  my $size = length($pdf_text);
  print "Content-Type: application/pdf\n";
  print "Content-Length: $size\n";
  print "$method\n";
  print "Content-Transfer-Encoding: binary\n\n"; # blank line to separate headers
  print $pdf_text;
}

sub invade_exteins{
  #takes a set of paired intein and extein trees as input
  #invades the extein trees with the inteins using a monte carlo chain
  #returns a nested hash for which the key is the extein tree, and has
  #data for which tips are invaded, and with which intein as nested data
  my $window = shift;
  my %tree_samples = %{my $hashref = shift};
  my %trees_with_invasion_data;
  foreach my $extein_tree (keys %tree_samples){
    my($extein_tree_object,@extein_tips)=readin_newick($extein_tree);
    my($intein_tree_object,@intein_tips)=readin_newick($tree_samples{$extein_tree});

    #get pairwise distance between all nodes in intein and extein tree
    my %extein_pairwise_distances= pairwise_patristic_distance($extein_tree_object,@extein_tips);
    my %intein_pairwise_distances= pairwise_patristic_distance($intein_tree_object,@intein_tips);

    #select a random intein and random extein to invade w/ said intein
    #then remove them from the list of valid invaders/invasion sites
    my $intein_tip = @intein_tips[rand($#intein_tips)];
    @intein_tips = remove_array_element($intein_tip,@intein_tips);
    my $extein_tip = @extein_tips[rand($#extein_tips)];
    @extein_tips = remove_array_element($extein_tip,@extein_tips);
    my %paired_sequence_archive;
    $paired_sequence_archive{$extein_tip}=$intein_tip;

    #take the random intein/extein pair as the starting point and begin the MC chain
    #effectively treats each infected tip as a seperate chain until all inteins have infected an extein
    until(!@intein_tips){
      my %paired_sequences = %paired_sequence_archive;
      foreach my $infected (keys %paired_sequences){
        #stop infecting if all inteins have infected
        next if(!@intein_tips);
        #get next intein tip to invade with
        my $closest_intein = get_smallest_distance($paired_sequences{$infected},\%intein_pairwise_distances);
        #gets all distances from current tip
        my %distance_from_current_extein = get_distances_from_tip($infected,\%extein_pairwise_distances);
        #returns accepted move/infection step
        $extein_tip = monte_carlo_chain($window,\%distance_from_current_extein,@extein_tips);
        #removes the intein that just invaded from list of inteins that need to infect a tip
        @intein_tips = remove_array_element($closest_intein,@intein_tips);
        #removes the just invaded extein as a valid move
        @extein_tips = remove_array_element($extein_tip,@extein_tips);
        $paired_sequence_archive{$extein_tip}=$intein_tip;
      }
    }

    #in theory, the %paired_sequence_archive hash now has the paired inteins and exteins\
    #now nest that data into the hash to be returned
    foreach my $extein (keys %paired_sequence_archive){
      $trees_with_invasion_data{$extein_tree}{"tips"}{$extein}=$paired_sequence_archive{$extein};
    }
    $trees_with_invasion_data{$extein_tree}{"intein_tree"}=$tree_samples{$extein_tree};
    #now loop to the top for all other trees in the set!
  }
  return(%trees_with_invasion_data);
}

sub remove_array_element{
  #takes value in array and array as inputs
  #returns array without that value
  my $index = 0;
  my $value_to_remove = shift;
  my @array = @_;
  $index++ until $array[$index] eq $value_to_remove;
  splice(@array, $index, 1);
  return(@array);
}


sub monte_carlo_chain{
  #takes a window size and hash of distances as inputs
  #returns a new tip that has been infected with an intein once the chain returns non-zero
  my $window_input = shift;
  my %hash_of_distances = %{my $hashref = shift};
  my @invaded_targets = @_;
  my $return_value = monte_carlo_step($window_input,\%hash_of_distances,@invaded_targets);
  if($return_value == 0){
    $return_value = monte_carlo_step($window_input,\%hash_of_distances,@invaded_targets);
  }
  else{
    return($return_value);
  }
}

sub get_distances_from_tip{
  #takes a nested hash of distances and key1
  #the hash is formatted as key1->key2->value
  #returns a hash of key2's with the associated distance to key1
  my $key1 = shift;
  my %distance_hash = %{my $hashref = shift};
  my %return_hash;
  foreach my $key (keys %distance_hash){
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
  my %distance_hash = %{my $hashref = shift};
  my $key_holder;
  foreach my $key (keys %distance_hash){
    next if($key ne $key1);
    my $smallest = "toggle";
    foreach my $key2 (keys %{$distance_hash{$key1}}){
      if(!$smallest eq "toggle"){
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
      $distances{$tip}{$tip_to_compare}=$node_distance;
    }
  }
  return(%distances);
}

sub monte_carlo_step{
  #takes a window width
  #and an array with associated distances from the current position
  #returns the new location
  my $window_width = shift;
  my %distance_of_tips = %{my $hashref = shift};
  my @invaded_exteins = @_;
  my $current_position = 0;

  #uniform random deviant
  my $uniform_random = rand();

  #creates a integer value based on the window size then gets new position proposal
  my $new_position = $current_position + int($uniform_random*$window_width);

  my $counter = 0;
  foreach my $key (sort {$distance_of_tips{$a} <=> distance_of_tips{$b}} keys %distance_of_tips){
    if($counter == $new_position){
      #first, check if the target has been invaded. If so, return 0
      if ( grep( /^$key$/, @invaded_exteins) ) {
        return(0);
      }
      #accept of reject the move based on an exponential distribution
      #using the distance from the old tip to the new tip as the X for the distribution
      my $acceptor = accept_or_deny_move($distance_of_tips{$key});
      if($acceptor == 1){
        return($key);
      }
      else{
        #return 0 if move not accepted
        return(0);
      }
    }
    else{
      $counter++;
      next;
    }
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

  #second uniform random
  my $uniform_random_deviant = rand();

  #compare uniform random to prob of acceptance. If lower, new position is accpeted!
  if($uniform_random_deviant <= $probability_of_acceptance){
    return(1);
  }
  else{
    return(0);
  }
}

sub simulate_sequences_from_tree{
  #takes an array of tree files as input
  #returns location of simulated sequence file in fasta format
  my $seq_type = shift;
  my $number_of_seqs = shift;
  my @tree_files = @_;
  my %trees;
  foreach my $file (@tree_files){
    my $tree_text = "";
    open(IN, "< $file");
    while(<IN>){
      $tree_text .= $_;
    }
    close IN;
    $trees{$file}=$tree_text;
  }

  #get user parameters for sequence sim
  print "Please provide inputs for $seq_type sequence simulation.\n
  Number of nucleotides per sequence?\n";
  my $num_of_nucs = <STDIN>;
  chomp $num_of_nucs;
  print "Alpha parameter for gamma distribution?\n";
  my $alpha = <STDIN>;
  chomp $alpha;
  print "How many gamma categories for the gamma distribution?\n";
  my $gamma_cats = <STDIN>;
  chomp $gamma_cats;

  #directory prep
  mkdir("simulated_sequences");
  chdir("simulated_sequences");
  my %return_hash;
  my $counter = 1;
  foreach my $tree (keys %trees){
    my $newick = $trees{$tree};
    #set up subdirectory to store results
    mkdir("$counter");

    #create random seed for simulation
    my $random_seed = int(rand(100000));
    if(0 == $random_seed % 2){
      $random_seed++;
    }

    #create simulation file
    #in the following formate
    #paml_format (uses 0 as default)
    #random_seed
    #number_of_seqs number_of_nucleotides number_of_replicates
    #tree_length (use -1 as the tree has absolute branch lengths)
    #tree_in_newick_format
    #alpha categories_for_gamma (gamma distribution)
    #model model_location (since using 2 -> WAG, need to provide location)
    #base frequencies table
    open(EVO, "+> $counter/MCaa.dat");
    print EVO
    "0\n
    $random_seed\n
    \n
    $number_of_seqs $num_of_nucs 1\n
    \n
    -1\n
    \n
    $newick\;\n
    \n
    $alpha $gamma_cats\n
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

    #now convert this file to fasta format
    #this might need correction......
    my $paml_file = "mc.paml";
    my $in  = Bio::AlignIO->new(-file   => $paml_file , -format => 'phylip', -interleaved => 1);
    my $out = Bio::AlignIO->new(-file   => ">$tree.fasta" , -format => 'fasta');
    while ( my $aln = $in->next_aln() ) {
      $out->write_aln($aln);
    }
    chdir("..");
    $return_hash{$tree}="simulated_sequences/$counter/$tree.fasta";
    $counter++;
  }
  chdir("..");
  return($num_of_nucs,%return_hash);
}

sub insert_sequences{
  #takes sequence and tree files for inteins and exteins
  #and the hash containing nested invasion data for each extein simulation
  #This sub then takes these inputs, and inserts the inteins into the exteins, and returns the files
  my $extein_length = shift;
  my %invasion_data = %{my $hashref1 = shift};
  my %extein_files = %{my $hashref2 = shift};
  my %intein_files = %{my $hashref3 = shift};

  #foreach extein tree file insert intein sequences into extein sequences in accordance with the invasion
  foreach my $invaded_extein_tree (keys %invasion_data){
    #select a random position to insert the intein into within the extein.
    #in order to be approximately similar to real life invasions, the invasion site
    #is placed within the 2nd or 3rd quartile of the sequence length
    my $quartile = int($extein_length/4);
    my $random_half = int(rand(($extein_length/2)));
    my $insertion_position = $random_half+$quartile;

    #get paire file information
    my $paired_intein_tree_file = $invasion_data{$invaded_extein_tree}{"intein_tree"};
    my $extein_sequence_file = $extein_files{$invaded_extein_tree};
    my $intein_sequence_file = $intein_files{$paired_intein_tree_file};

    #readin fasta files to memory
    my %extein_sequences = readin_fasta($extein_sequence_file);
    my %intein_sequences = readin_fasta($intein_sequence_file );

    #foreach extein sequence tip associated with the tree, insert the associated intein sequence
    foreach my $invaded_extein_tip (keys %{$invasion_data{$invaded_extein_tree}{"tips"}}){
      #variable name makes sense shortly
      my $extein_left_sequence = $extein_sequences{$invaded_extein_tip};
      my $extein_right_sequence = substr $extein_left_sequence, 0, $insertion_position, '';
      my $intein_name = $invasion_data{$invaded_extein_tree}{"tips"}{$invaded_extein_tip};
      my $intein_sequnece = $intein_sequences{$intein_name};
      my $extein_w_intein_sequence = $extein_left_sequence.$intein_sequnece.$extein_right_sequence;
      delete($extein_sequences{$invaded_extein_tip});
      my $new_tip_name = $invaded_extein_tip."_w_intein";
      $extein_sequences{$new_tip_name}=$extein_w_intein_sequence;
    }

    #print invaded and uninvaded extein seuqneces to file
    mkdir("invaded_sequences");
    open(OUT, "+> $extein_sequence_file.invaded");
    foreach my $accession (keys %extein_sequences){
      print OUT "\>$accession\n";
      print OUT "$extein_sequences{$accession}\n";
    }
    close OUT;
  }
}

sub readin_fasta{
  my $infile = shift;
  my $accession="";
  my %sequences;
  open(IN, "< $infile");
  while(<IN>){
    chomp;
    if($_=~/\>/){
      $accession=$_;
      $sequences{$accession}="";
    }
    else{
      $sequences{$accession}.=$_;
    }
  }
  close IN;
  return(%sequences);
}
