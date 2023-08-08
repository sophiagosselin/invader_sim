#!/usr/bin/perl -w
use strict;
use warnings;

#this will get me some statistics

my $directory = "ice_blast_runs";

open(my $stats, "+> indv_sample_results.txt");
#header
print $stats "Sample_Number\tIdentity_Value\tE_Value\tMatches_Found\tInteins_Found\tNon_Intein_Matches\tAverage_Length_of_Intein_Found\n";

MAIN();

close $stats;

sub MAIN{
  opendir DIR, $directory ;
  while(my $dir = readdir DIR){
    # $dir will be the sample number
    if(-d "$directory\/$dir"){
      opendir DIR_SUB, "$directory\/$dir";
      while (my $dir_sub = readdir DIR_SUB){
        # $dir_sub will be the parameter directory
        if(-d "$directory\/$dir\/$dir_sub"){
          opendir DIR_SUB_SUB, "$directory\/$dir\/$dir_sub";
          while (my $dir_sub_sub = readdir DIR_SUB_SUB){
            # $dir_sub_sub is the directory where ice blast ran
            if($dir_sub_sub eq "output"){
              #print "File is $directory\/$dir\/$dir_sub\/$dir_sub_sub\/all_matches.fasta\n";
              my %output_sequences=READIN_FASTA("$directory\/$dir\/$dir_sub\/$dir_sub_sub\/all_matches.fasta");
              my $seq_avg = my $intein_counter = my $match_counter = my $fp_matches =0;
              foreach my $asc (keys %output_sequences){
                $match_counter++;
                if($asc=~/w_intein/){
                  $intein_counter++;
                  $seq_avg+=length($output_sequences{$asc});
                }
                else{
                  $fp_matches++;
                }
              }
              $seq_avg=$seq_avg/$intein_counter;
              my($identity,$eval)=($dir_sub=~/(.*?)(1e.*)/);
              #print "Troubleshooting\nSample number is $dir\nIdentity is $identity\nEval is $eval\nNum of Matches is $match_counter\nNum of inteins is $intein_counter\nNum of FPs is $fp_matches\n";
              print $stats "$dir\t$identity\t$eval\t$match_counter\t$intein_counter\t$fp_matches\t$seq_avg\n";
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
    }
    else{
      next;
    }
  }
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
