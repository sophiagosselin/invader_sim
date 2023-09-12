#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;
use Getopt::Long;
use Bio::TreeIO;
use Bio::AlignIO;
use Bio::Tree::Node;
use Bio::Tree::TreeFunctionsI;
use Cwd;
use threads;
use threads::shared;
use Data::Dumper;


$| = 1;
#future coding notes:
#might be worth adding a scratch file to solve the hanging issue
#should rewrite

#GLOBALS
#naiive mode boolean: 0=naiive mode, 1=startup file inputs
my $naiive_mode = 0;
my (%extein_tree_params, %intein_tree_params, %sample_params, %extein_seq_params, %intein_seq_params);
my %sample_getopt = ("sn"=>"sn=i","ws"=>"ws=i","evo"=>"evo=s","sub"=>"sub:i");
my %tree_sim_getopt = ("sp"=>"sp=i","b"=>"b=s","d"=>"d=s","sf"=>"sf=s","m"=>"m=s");
my %sequence_sim_getopt = ("nn"=>"nn=s","a"=>"a=s","cat"=>"cat=s","m"=>"m=s");

#Code Start

STARTUP();
MAIN();

sub STARTUP{
  #check for naive mode vs. file with parameter inputs.
  #also checks for help call
  #if naive mode is called, then proceed. If not, load parameters from file
  #parameter file should be in the following format

  #first, check for help call
  my $help = "0";
  GetOptions('help+' => \$help, 'h+' => \$help);
  if($help ne "0"){
    die
"\n\n\n****************************************************************************

INVADER-SIM (Intein Invasion Sequence Simulator) v1.1.0\n

Author: Gosselin Sophia
Bug reporting: https://github.com/sophiagosselin/evolver_pipeline

This program will simulate extein and intein sequences according to user parameters.
Then using a monte carlo based approach, invade the extein sequences with the intein sequences.
Returns simulated samples of extein nucleotide sequences invaded by intein sequences.

Usage: perl invader_sim.pl


The program can be run in naiive or prepared mode.
If the user wishes to be guided through the simulation process with prompts, then simply launch the program as above.


Else, a parameter file (is.param) can be provided using the following parameters:

\#These parameters apply to phylogeny simulations
-sp -> number of species.
-b -> birth rate.
-d -> death rate.
-sf -> sample fraction.
-m -> mutation rate.\n

\#These parameters apply to sequence simulations
-nn -> number of nucleotides.
-a -> alpha parameter for gamma distribution
-cat -> number of rate categories for the gamma distribution
-m -> model number: (0: Poisson, 1: Proportional)

\#These parameters are general to the code
-sn -> number of simulations to create
-ws -> window size for MC chain (how many tips away can an intein jump)
-evo -> string to call paml-evolver (likely paml-evolver, or evolver)
-sub -> (OPTIONAL) Number of intein sequences to invade with out of the set of simulated sequences. Useful if you want to simulate X inteins but only invade with Y of them.

Note that the is.param must be formated as:
extein phylogeny simulation parameters
intein phylogeny simulation parameters
global params
extein sequence simulation parameters
intein sequence simulation parameters

An example is.param file:
-sp 50 -b 1 -d 1 -sf .1 -m .01
-sp 10 -b 1 -d 1 -sf .1 -m .1
-sn 100 -ws 10 -evo evolver -sub 5
-nn 100 -a .5 -cat 4 -m 0
-nn 50 -a 1 -cat 4 -m 1 \n\n";
  }
  if(-e "is.param"){
    $naiive_mode = 1;
    my $counter = 0;
    open(IN, "< is.param");
    while(<IN>){
      chomp;
      next if($_!~/\-/);
      if($counter != 5){
          if($counter == 0){
            %extein_tree_params = parse_and_check_inputs(1,"",\%tree_sim_getopt,$_);
          }
          elsif($counter == 1){
            %intein_tree_params = parse_and_check_inputs(1,"",\%tree_sim_getopt,$_);
          }
          elsif($counter ==2){
            %sample_params = parse_and_check_inputs(1,"",\%sample_getopt,$_);
          }
          elsif($counter == 3){
            %extein_seq_params = parse_and_check_inputs(1,"",\%sequence_sim_getopt,$_);
          }
          elsif($counter == 4){
            %intein_seq_params = parse_and_check_inputs(1,"",\%sequence_sim_getopt,$_);
          }
          $counter++;
      }
      else{
        last;
      }
    }
  }
  else{
    $sample_params{"evo"} = parse_and_check_inputs(0,"\n\n****************************************************************************\nPlease input the call for paml-evolver on this machine.\nThis is likely paml-evolver, or evolver.\n\n","",0);
  }
}

sub MAIN{
  #get user inputs and test until they are satisfied
  if($naiive_mode == 0){
    %extein_tree_params = test_parameters("extein");
    %intein_tree_params = test_parameters("intein");

    #once the user is satisfied with all parameters, ask for number of experiments to run
    $sample_params{"sn"} = parse_and_check_inputs(0,"\n\n****************************************************************************\nParameters have been set; ready to run experiment.\n How many samples would you like to create?\n\n","",0);
  }

  #create the number of trees for inteins and exteins based on the number of experiments requested
  my @extein_phylogenies = simulate_n_trees("extein",$sample_params{"sn"},\%extein_tree_params);
  my @intein_phylogenies = simulate_n_trees("intein",$sample_params{"sn"},\%intein_tree_params);

  #important to make copy for array eating process
  my @intein_phylogenies_for_pairing = @intein_phylogenies;
  #pair each extein tree with a random intein tree. Each pairing constitutes a sample
  my %paired_intein_extein_trees;
  foreach my $extein_tree (@extein_phylogenies){
    my $random_index = rand(@intein_phylogenies_for_pairing);
    my $random_intein_tree = splice(@intein_phylogenies_for_pairing,$random_index,1);
    $paired_intein_extein_trees{$extein_tree}=$random_intein_tree;
  }

  if($naiive_mode == 0){
    #ask the user for max jump size when invading from one extein tip to another
    $sample_params{"ws"} = parse_and_check_inputs(0,"\n\n****************************************************************************\nPreparing intein invasion simulation.\nWhat window size would you like to constrain intein invasions to?\n\nFor a simulation with n exteins, and x inteins, this window (w) should follow this formula:\nw <= n-x\n\nFor example, in a tree with 10 tips, a window of 5 would allow the intein to invade any of the closest 5 tips to its current position.\n\nNote that exceeding the size constraint of the above equation will lead to unconstrained invasion (as, eventually, all tips will be available for invasion)\n\n","",0);
  }

  #simulates intein invasion
  #returns a nested hash in the following format:
  #$extein tree -> {"tips"/"intein_tree"} -> {$extein_tips/$intein_file_name} -> $invading_intein
  my %paired_invasion_data = invade_exteins($sample_params{"ws"},\%paired_intein_extein_trees);

  #now simulate sequence evolution based on the phylogenies
  #returns a hash of the sequence file associated with the tree file
  my ($extein_size,%paired_extein_tree_and_sequence_files) = simulate_sequences_from_tree("extein",$extein_tree_params{"sp"},\%extein_seq_params,@extein_phylogenies);
  my ($intein_size,%paired_intein_tree_and_sequence_files) = simulate_sequences_from_tree("intein",$intein_tree_params{"sp"},\%intein_seq_params,@intein_phylogenies);

  #now invade extein sequences with intein sequences
  insert_sequences($extein_size,\%paired_invasion_data,\%paired_extein_tree_and_sequence_files,\%paired_intein_tree_and_sequence_files);

  print "\n\n****************************************************************************\n\nSimulation finished! Find the fasta files of invaded extein sequences with their un-invaded counterparts in the 'invaded_sequences' directory.\n\n";
}

sub test_parameters{
  #takes tree type as input, returns user determined inputs for that tree type.
  my $tree_type = shift;
  my $user_check = "N";
  my $input_primer = "\n\n****************************************************************************
Please input parameters for the $tree_type phylogeny now!
Remember that the extein tree should have at least 2x the amount of taxa as the intein tree.
(This is critical for the method in which this software simulates intein invasion)\n
Parameters should be in the following format:\n
sp - number of species.
b - birth rate.
d - death rate.
sf - sample fraction.
m - mutation rate.\n
Example: -sp 50 -b 1 -d 1 -sf .1 -m .01\n\n";
  my %test_parameters = parse_and_check_inputs(1,$input_primer,\%tree_sim_getopt,0);
  #begin parameter testing for phylogeny of interest
  print "\n\nBeginning extein simulation parameters test.\n Testing will conclude once the user is happy with the tree and set parameters.\n\n";
  %test_parameters = parameter_testing_loop($user_check,$tree_type,$input_primer,\%test_parameters);
  #save parameters settled on
  return(%test_parameters);
}

sub parse_and_check_inputs{
  #takes a boolean to determine whether to parse a single or multiple inputs
  #0 for single with no flags, 1 for multiple inputs with flags
  #takes a string as the message for the user and primes them to input options
  #takes a hash of keys (which is the input tag the user is inputting)
  #Example: %hash{"input_name"}="parse_string_for_GetOptions" <- see GetOptions for info on what each type is
  #then splits apart that response, and assigns values to global variables
  #additionally can take a string as input instead of the user provided input
  my $string_or_hash = shift;
  my $message = shift;
  my $hashref = shift;
  my %input_types;
  if(defined $hashref && $hashref ne ''){
    %input_types = %{$hashref};
  }
  my $new_options = shift;

  #prime user for input if necessary
  if($new_options eq "0"){
    print($message);
    $new_options = <STDIN>;
    chomp $new_options;
  }

  #return string if that is all that is requested
  if($string_or_hash == 0){
    return($new_options);
  }
  elsif($string_or_hash == 1){
    # Split the user-provided options into individual arguments
    my @options = split /\s+/, $new_options;

    # Add the options as separate elements in @ARGV
    @ARGV = @options;

    #parse options
    my %parsed_inputs;
    #this syntax might not work...
    my @opts;
    foreach my $input_key (keys %input_types) {
      push @opts, $input_types{$input_key} => \$parsed_inputs{$input_key};
    }
    GetOptions(@opts);

    #check input variables
    print "\n\nChecking Parameters\n\n";
    foreach my $input (keys %parsed_inputs){
      if(!defined $parsed_inputs{$input} || $parsed_inputs{$input} eq 'undef'){
        if($naiive_mode == 0){
          print "$input is not defined. Please define it now:\n";
          my $new_parameter = <STDIN>;
          chomp $new_parameter;
          $parsed_inputs{$input}=$new_parameter;
        }
        elsif($input eq "sub"){
          $parsed_inputs{$input}=0;
          next;
        }
        else{
          die "$input not defined\n";
        }
      }
      else{}
    }
    return(%parsed_inputs);
  }
  else {
    die "Invalid value for \$string_or_hash: $string_or_hash";
  }
}

sub parameter_testing_loop{
  #takes a y/n toggle, tree extsension and prossibly tree object as inputs
  #loops through tree building process until the user is satisfied based on the
  #global parameters
  #returns nothing, but the final chosen parameters will be in the global hash %inputs
  my $toggle = shift;
  my $tree_extension = shift;
  my $input_msg = shift;
  my %params = %{my $hashref = shift};
  if(uc($toggle) eq "N"){
    my $evolver_tree = simulate_1_tree(\%params);
    my $pdf_tree = write_tree($evolver_tree);
    $toggle = parse_and_check_inputs(0,"\n\n****************************************************************************\nDoes this simulated $tree_extension tree look good to you [Y/N]?\n","",0);
    if(uc($toggle) ne "Y"){
      %params = parse_and_check_inputs(1,$input_msg,\%tree_sim_getopt,0);
      %params = parameter_testing_loop($toggle,$tree_extension,$input_msg,\%params);
    }
    elsif(uc($toggle) eq "Y"){
      return(%params);
    }
  }
  elsif(uc($toggle) eq "Y"){
    return(%params);
  }
  return(%params);
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

  my @loop = (1..$simulation_num);
  foreach my $num (@loop){
    my $phylogeney = simulate_1_tree(\%inputs_for_tree);
    my $path_to_file = "$file_handle\/$file_handle\_$num\.tree";
    copy($phylogeney,$path_to_file);
    push(@tree_files,$path_to_file);
    unlink $phylogeney;
  }
  return(@tree_files);
}

sub simulate_1_tree{
  #simulated a tree based on parameters provided
  #first generates a random seed and ensures it is an odd value
  my %tree_parameters = %{my $hashref = shift};

  my $random_tree_seed = int rand(1000000000);
  if(0 == $random_tree_seed % 2){
    $random_tree_seed++;
  }

  #get input params from %inputs
  my $species_number=$tree_parameters{'sp'};
  my $tree_number=1;
  my $birth_rate=$tree_parameters{'b'};
  my $death_rate=$tree_parameters{'d'};
  my $sample_fraction=$tree_parameters{'sf'};
  my $mutation_rate=$tree_parameters{'m'};

  open(my $evolver_pipe, '|-', $sample_params{"evo"}) or die "Could not open pipe to paml-evolver: $!\n";
  print $evolver_pipe "1\n"; #uses random unrooted tree option
  print $evolver_pipe "$species_number\n"; #species number in tree
  print $evolver_pipe "1 $random_tree_seed\n"; #number of trees and random seed
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

  #Note: this should be rewritten at some point such that the node objects are created
  #and then ID's and distances between ID's are extracted. Such that the MC chain
  #does not interact with the objects but just their properties.

  #takes a set of paired intein and extein trees as input
  #invades the extein trees with the inteins using a monte carlo chain
  #returns a nested hash for which the key is the extein tree, and has
  #data for which tips are invaded, and with which intein as nested data
  my $window = shift;
  my %tree_samples = %{my $hashref = shift};
  my %trees_with_invasion_data;
  foreach my $extein_tree (keys %tree_samples){
    #readin newick, get tree object, and node objects
    print "\n\nReading in trees to memory for paired sample $extein_tree and $tree_samples{$extein_tree}.\n\n";
    my($extein_tree_object,$array_ref)=readin_newick($extein_tree);
    my @extein_tip_objects = @{$array_ref};
    my($intein_tree_object,$array_ref2)=readin_newick($tree_samples{$extein_tree});
    my @intein_tip_objects = @{$array_ref2};

    #get pairwise distance between all nodes in intein and extein tree
    #importantly this hash is structured %hash{"nodeID"}{"nodeID"}=distance_value
    print "Calculating pairwise patristic distance between tips.\n\n";
    my ($array_of_hashes_ref)=pairwise_distance_high($extein_tree_object,@extein_tip_objects);
    my @array_of_hashes = @{$array_of_hashes_ref};
    my (%extein_pairwise_distances) = MERGE_NESTED_HASHES(@array_of_hashes);

    my ($array_of_hashes_ref2)=pairwise_distance_high($intein_tree_object,@intein_tip_objects);
    my @array_of_hashes2 = @{$array_of_hashes_ref2};
    my (%intein_pairwise_distances) = MERGE_NESTED_HASHES(@array_of_hashes2);

    #clear memory
    $array_of_hashes_ref="";
    $array_of_hashes_ref2="";
    @array_of_hashes=();
    @array_of_hashes2=();

    #gets arrays of tip ids for inteins and exteins
    my @intein_tip_ids;
    foreach my $node_obj (@intein_tip_objects){
      my $id_of_node = get_id($node_obj,@intein_tip_objects);
      push(@intein_tip_ids,$id_of_node);
    }
    my @extein_tip_ids;
    foreach my $node_obj (@extein_tip_objects){
      my $id_of_node = get_id($node_obj,@extein_tip_objects);
      push(@extein_tip_ids,$id_of_node);
    }

    #if subset mode is enabled, randomly select the subset of intein tips now.
    my @invaded_inteins;
    if($sample_params{"sub"}!=0){
      my @random_subsample;
      my @intein_tips_for_sub=@intein_tip_ids;
      for(my $sub_counter=0; $sub_counter<$sample_params{"sub"}; $sub_counter++){
        my $rand_tip = splice(@intein_tips_for_sub, rand @intein_tips_for_sub, 1);
        push(@random_subsample,$rand_tip);
      }
      #pushes unused intein to the invaded list, such that they cannot be used again
      @invaded_inteins=@intein_tips_for_sub;
      @intein_tip_ids=@random_subsample;
    }
    else{}

    #select a random intein and random extein to invade w/ said intein
    #then remove them from the list of valid invaders/invasion sites
    #then save pairing of object refs
    print "Preparing starting point for invasion.\n\n";
    my (@invaded_exteins,%paired_sequence_archive);
    my $intein_tip = splice(@intein_tip_ids, rand @intein_tip_ids, 1);
    my $extein_tip = splice(@extein_tip_ids, rand @extein_tip_ids, 1);
    my %paired_extein_intein_tip_ids = ($extein_tip => $intein_tip);
    $paired_sequence_archive{$extein_tip}=$intein_tip;
    push(@invaded_exteins,$extein_tip);
    push(@invaded_inteins,$intein_tip);
    #print "First pair is Extein: $extein_tip Intein: $intein_tip\n";

    #take the random intein/extein pair as the starting point and begin the MC chain
    #effectively treats each infected tip as a seperate chain until all inteins have infected an extein
    print "Starting Monte-Carlo chain based invasion simulation.\n\n";
    until(!@intein_tip_ids){
      my %itterating_paired_tips = %paired_extein_intein_tip_ids;
      foreach my $extein_tip_for_work (keys %itterating_paired_tips){
        #print "Extein to invade from: $extein_tip_for_work\n";
        #stop infecting if all inteins have infected
        next if(!@intein_tip_ids);

        #get next intein tip to invade with.
        #inputs ref to intein that invaded current extein of foreach loop.
        my $new_intein_id = get_smallest_distance($itterating_paired_tips{$extein_tip_for_work},\%intein_pairwise_distances,@invaded_inteins);
        #print "Closest intein: $new_intein_id\n";

        #gets all distances from current tip excluding those previously searched
        my %distance_from_current_extein = get_distances_from_tip($extein_tip_for_work,\%extein_pairwise_distances,@invaded_exteins);

        #returns accepted move/infection step
        my $new_extein_id= monte_carlo_chain($window,\%distance_from_current_extein);
        #print "Extein to invade into $new_extein_object\n";
        #defref object

        #retrieve original object, and use it to get ID's for downstream use
        #print "$new_extein_id\t$new_intein_id\n";
        $paired_sequence_archive{$new_extein_id}=$new_intein_id;

        #removes the intein that just invaded from list of inteins that need to infect a tip
        @intein_tip_ids = remove_array_element($new_intein_id,@intein_tip_ids);
        #removes the just invaded extein as a valid move
        @extein_tip_ids = remove_array_element($new_extein_id,@extein_tip_ids);

        push(@invaded_exteins,$new_extein_id);
        push(@invaded_inteins,$new_intein_id);

        $paired_extein_intein_tip_ids{$new_extein_id}=$new_intein_id;
      }
    }

    print "Invasion for paired sample complete.\n\n";
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

sub MERGE_NESTED_HASHES{
	my @hashes_to_merge = @_;
	my %new_hash;
	foreach my $hashref (@hashes_to_merge){
		my %temp_hash = %{$hashref};
		%new_hash = (%new_hash,%temp_hash);
	}
	my $test = (\%new_hash);
	return(%new_hash);
}

sub get_id{
  #takes object array and ref to node object
  #returns tip id of node object
  my $obj_ref = shift;
  my @obj_array = @_;
  my $index = 0;
  foreach my $obj (@obj_array){
    #print "Array includes $value\n";
    if($obj_array[$index] eq $obj_ref){
      my $id = $obj->Bio::Tree::Node::id();
      return($id);
    }
    else{
      $index++;
    }
  }
}

sub remove_array_element{
  #takes value in array and array as inputs
  #returns array without that value
  my $value_to_remove = shift;
  my @array = @_;
  my $index = 0;
  my $toggle =0;
  foreach my $value (@array){
    #print "Array includes $value\n";
    if($array[$index] eq $value_to_remove){
      splice(@array,$index,1);
      $toggle = 1;
    }
    else{
      $index++;
    }
  }
  if($toggle == 0){
    die "\n\nCould not find $value_to_remove\n\n";
  }
  return(@array);
}

sub get_distances_from_tip{
  #takes a nested hash of distances and key1
  #the hash is formatted as key1->key2->value
  #can also take an array of values to skip if requested
  #returns a hash of key2's with the associated distance to key1
  my $key1 = shift;
  my %distance_hash = %{my $hashref = shift};
  my %return_hash;
  my @keys_to_skip = @_;
  foreach my $key (keys %distance_hash){
    foreach my $key2 (keys %{$distance_hash{$key1}}){
      next if($key2 eq $key);
      my $toggle = 0;
      foreach my $test (@keys_to_skip){
        if($key2 eq $test){
          $toggle =1;
        }
      }
      next if($toggle == 1);
      $return_hash{$key2}=$distance_hash{$key1}{$key2};
    }
  }
  return(%return_hash);
}

sub get_smallest_distance {
  #takes a nested hash of distances and key1
  #the hash is formatted as key1->key2->value
  #can also take an array of values to ignore
  #returns the key2 with the smallest distance to key1
    my $key_to_find = shift;
    my %distance_hash = %{ shift() };
    my @keys_to_skip = @_;
    my $key_holder;
    my $smallest;
    foreach my $key1 (keys %distance_hash) {
      #print "Key1 $key1\n";
        next if ($key1 ne $key_to_find);

        foreach my $key2 (keys %{ $distance_hash{$key1} }) {
          #print "Key2 $key2\n";
            next if ($key2 eq $key_to_find);
            my $toggle = 0;
            foreach my $test (@keys_to_skip){
              if($test eq $key2){
                $toggle = 1;
              }
            }
            next if($toggle == 1);
            if (!defined $smallest || $distance_hash{$key1}{$key2} < $smallest) {
                $smallest = $distance_hash{$key1}{$key2};
                $key_holder = $key2;
            }
        }
    }
    #print "I'm returning $key_holder\n";
    return $key_holder;
}

sub readin_newick{
  #readin newick format tree file
  #return a tree object and an array of tips
  my $tree_file = shift;
  my $tree_input = Bio::TreeIO->new(-file => $tree_file,
                            -format => "newick");
  my $tree_object = $tree_input->next_tree;
  my @tree_taxa = $tree_object->get_leaf_nodes;
  return($tree_object,\@tree_taxa);
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

sub pairwise_distance_high{
  #takes an array of tree tip nodes, and a tree object as inputs
  #returns a nested hash of object reference keys paired to their distance from the node of interest
  my $tree_object = shift;
  my @array_of_tips = @_;
  my (@threads,@return_values);
  my $subroutine="pairwise_distance_low";
  foreach my $tip (@array_of_tips){
    my($thread)=threads->create(\&$subroutine,$tree_object,$tip,@array_of_tips);
    push(@threads,$thread);
  }
  foreach (@threads){
    my $return_value = $_->join();
    push(@return_values,$return_value);
  }
  #print "Finished pairwise sub\n";
  return(\@return_values);
}

sub pairwise_distance_low{
  #takes an array of tree tip nodes, and a tree object as inputs
  #returns a nested hash of object reference keys paired to their distance from the node of interest
  my $tree_object = shift;
  my $tip = shift;
  my %distances;
  my @array_of_tips = @_;
  foreach my $tip_to_compare (@array_of_tips){
    my $node_distance = get_patristic_distance($tip,$tip_to_compare,$tree_object);
    my $tipid1=get_id($tip,@array_of_tips);
    my $tipid2=get_id($tip_to_compare,@array_of_tips);
    $distances{$tipid1}{$tipid2}=$node_distance;
  }
  return(\%distances);
}

sub monte_carlo_chain{
  #takes a window size and hash of distances as inputs
  #returns a new tip that has been infected with an intein once the chain returns non-zero
  my $window_input = shift;
  my %hash_of_distances = %{my $hashref = shift};
  my $return_value = monte_carlo_step($window_input,\%hash_of_distances);
  if($return_value eq "0"){
    $return_value = monte_carlo_chain($window_input,\%hash_of_distances);
  }
  else{
    return($return_value);
  }
  return($return_value);
}

sub monte_carlo_step{
  #takes a window width
  #and an array with associated distances from the current position
  #returns the new location
  my $window_width = shift;
  my %distance_of_tips = %{my $hashref = shift};

  #uniform random deviant
  my $uniform_random = rand();

  #creates a integer value based on the window size then gets new position proposal
  my $new_position = int($uniform_random*$window_width);

  my $counter = 0;
  foreach my $key (sort {$distance_of_tips{$a} <=> $distance_of_tips{$b}} keys %distance_of_tips){
    if($counter == $new_position){
      #accept of reject the move based on an exponential distribution
      #using the distance from the old tip to the new tip as the X for the distribution
      my $acceptor = accept_or_deny_move($distance_of_tips{$key});
      if($acceptor == 1){
        return($key);
      }
      else{
        #return 0 if move not accepted
        return("0");
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
  my %sequence_sim_params = %{my $hashref = shift};
  my @tree_files = @_;
  my $counter = 1;
  my $kappa = "";
  my (%return_hash);
  #directory prep
  if(!-d "simulated_sequences"){
    mkdir("simulated_sequences");
  }
  mkdir("simulated_sequences/$seq_type");

  #get user parameters for sequence sim
  if($naiive_mode == 0){
    %sequence_sim_params = parse_and_check_inputs(1,"\n\n****************************************************************************\n\nPlease provide inputs for $seq_type sequence simulation.\n
    The follwing inputs are required:
    -nn Number of nucleotides in simulated sequences
    -a alpha parameter for gamma distribution
    -cat number of rate categories for the gamma distribution
    -m model number. Options are: (0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV)
    -k kappa value. Only necessary if models 1 or 4 are selected.\n
    Example input: -nn 100 -a .5 -cat 4 -m 7\n\n",\%sequence_sim_getopt,0);
  }

  else{}

  #start simulating
  foreach my $file (@tree_files){
    #first, readin tree file
    my $newick = "";
    my $dir = getcwd;
    print "$file is the current simulation target\n";
    open(IN, "< $file") or die "$file could not be opened! Currently in $dir for WD\n";
    while(<IN>){
      $newick .= $_;
    }
    close IN;

    #check if newick has ; at the end. If not, add one
    chomp $newick;
    if($newick!~/.*\;/){
      $newick.=";";
    }

    #make and enter directory for simulation
    mkdir("simulated_sequences/$seq_type/$counter");
    chdir("simulated_sequences/$seq_type/$counter");

    #create random seed for simulation. Per PAML this must be an odd #
    my $random_seed = int(rand(100000));
    if(0 == $random_seed % 2){
      $random_seed++;
    }

    #note: Ignore the user manual for paml.
    #It incorrectly places the model before the gamma distribution
    #create simulation file in the following format:
    #paml_format (uses 0 as default)
    #random_seed
    #number_of_seqs number_of_nucleotides number_of_replicates
    #tree_length (use -1 as the tree has absolute branch lengths)
    #tree_in_newick_format
    #model # (0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV)
    #kappa/rate parameters in model
    #alpha + categories_for_gamma (gamma distribution)
    #base frequencies table
    open(EVO, "+> MCaa.dat");
    print EVO
    "0\n$random_seed\n$number_of_seqs $sequence_sim_params{'nn'} 1\n-1\n$newick\n\n$sequence_sim_params{'a'} $sequence_sim_params{'cat'}\n$sequence_sim_params{'m'}\n0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05\n0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05\n\nA R N D C Q E G H I\nL K M F P S T W Y V\n\n\/\/ end of file";
    close EVO;

    #call evolver and simulate. OUT is mc.paml in phylip format
    open(my $evolver_pipe, '|-', $sample_params{"evo"}) or die "Could not open pipe to paml-evolver: $!\n";
    print $evolver_pipe "7";
    close $evolver_pipe;

    #now convert this file to fasta format
    #this might need correction......
    my $paml_file = "mc.paml";
    my $in  = Bio::AlignIO->new(-file   => $paml_file, -format => 'phylip', -interleaved => 1);
    my $out = Bio::AlignIO->new(-file   => ">sample_$counter.fasta", -format => 'fasta');
    while ( my $aln = $in->next_aln() ) {
      $out->write_aln($aln);
    }

    #return location of sequence file
    $return_hash{$file}="simulated_sequences/$seq_type/$counter/sample_$counter.fasta";
    $counter++;
    chdir("..");
    chdir("..");
    chdir("..");
  }
  return($sequence_sim_params{"nn"},%return_hash);
}

sub insert_sequences{
  #takes sequence and tree files for inteins and exteins
  #and the hash containing nested invasion data for each extein simulation
  #This sub then takes these inputs, and inserts the inteins into the exteins, and returns the files
  my $extein_length = shift;
  my %invasion_data = %{my $hashref1 = shift};
  my %extein_files = %{my $hashref2 = shift};
  my %intein_files = %{my $hashref3 = shift};
  mkdir("invaded_sequences");
  print "Now inserting intein sequences into extein sequences.\n\n";

  #foreach extein tree file insert intein sequences into extein sequences in accordance with the invasion
  foreach my $invaded_extein_tree (keys %invasion_data){
    #select a random position to insert the intein into within the extein.
    #in order to be approximately similar to real life invasions, the invasion site
    #is placed within the 2nd or 3rd quartile of the sequence length
    my $quartile = int($extein_length/4);
    my $random_half = int(rand(($extein_length/2)));
    my $insertion_position = $random_half+$quartile;

    #get paired file information
    my $paired_intein_tree_file = $invasion_data{$invaded_extein_tree}{"intein_tree"};
    my $extein_sequence_file = $extein_files{$invaded_extein_tree};
    my $intein_sequence_file = $intein_files{$paired_intein_tree_file};
    #print "Paired information: $invaded_extein_tree\t$paired_intein_tree_file\t$extein_sequence_file\t$intein_sequence_file\n\n";

    #readin fasta files to memory
    my %extein_sequences = readin_fasta($extein_sequence_file);
    my %intein_sequences = readin_fasta($intein_sequence_file );

    #foreach extein sequence tip associated with the tree, insert the associated intein sequence
    foreach my $invaded_extein_tip (keys %{$invasion_data{$invaded_extein_tree}{"tips"}}){
      #print "Tip to invade $invaded_extein_tip\n";
      #variable name makes sense shortly
      my $extein_right_sequence = $extein_sequences{$invaded_extein_tip};
      my $extein_left_sequence = substr $extein_right_sequence, 0, $insertion_position, '';
      my $intein_name = $invasion_data{$invaded_extein_tree}{"tips"}{$invaded_extein_tip};
      my $intein_sequnece = $intein_sequences{$intein_name};
      my $extein_w_intein_sequence = $extein_left_sequence.$intein_sequnece.$extein_right_sequence;
      delete($extein_sequences{$invaded_extein_tip});
      my $new_tip_name = $invaded_extein_tip."_w_intein_".$intein_name;
      $extein_sequences{$new_tip_name}=$extein_w_intein_sequence;
    }

    #print invaded and uninvaded extein seuqneces to file
    my ($file_descriptor1)=($extein_sequence_file=~/.*\/(.*?)\.fasta/);
    my ($file_descriptor2)=($intein_sequence_file=~/.*\/(.*?)\.fasta/);
    open(OUT, "+> invaded_sequences/extein\_$file_descriptor1\_intein\_$file_descriptor2.fasta");
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
      #this regex may need change if paml changes.
      #Currently this does retrieve the correct tip values
      my ($no_carrot) = ($_=~/\>(.*?)\/.*/);
      $accession=$no_carrot;
      $sequences{$accession}="";
    }
    else{
      $sequences{$accession}.=$_;
    }
  }
  close IN;
  return(%sequences);
}
