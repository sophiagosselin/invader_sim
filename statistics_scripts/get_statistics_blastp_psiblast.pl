#!/usr/bin/perl -w
use strict;
use warnings;

#this will get me some statistics

my $directory = "ice_blast_runs";

open(my $blastp_stats, "+> blastp_sample_results.txt");
open(my $psiblast_stats, "+> psiblast_sample_results.txt");

#header
print $blastp_stats "Sample_Number\tE_Value\tMatches_Found\tInteins_Found\tNon_Intein_Matches\tAverage_Length_of_Intein_Found\n";
print $psiblast_stats "Sample_Number\tE_Value\tMatches_Found\tInteins_Found\tNon_Intein_Matches\tAverage_Length_of_Intein_Found\n";

MAIN();

close $blastp_stats;
close $psiblast_stats;

sub MAIN{
  opendir DIR, $directory ;
  while(my $sample_number = readdir DIR){
    if(-d "$directory\/$sample_number"){
      opendir DIR_SUB, "$directory\/$sample_number";
      while (my $sample_contents = readdir DIR_SUB){
        
        #for the blastp experiment
        if($sample_contents eq "blastp"){
          opendir BLASTP, "$directory\/$sample_number\/$sample_contents";
          while(my $blastp_runs = readdir BLASTP){
            opendir BLASTP_EVAL, "$directory\/$sample_number\/$sample_contents\/$blastp_runs";
            while(my $blastp_sample = readdir BLASTP_EVAL){
              
              #directory of interest
              if($blastp_sample eq "blastp.blast"){
                #gets statistics
                my($blastp_matches,$blastp_inteins,$blastp_avg)=PARSE_BLAST("$directory\/$sample_number\/$sample_contents\/$blastp_runs\/$blastp_sample");
                my $non_inteins = $blastp_matches-$blastp_inteins;

                #print to file
                print $blastp_stats "$sample_number\t$blastp_runs\t$blastp_matches\t$blastp_inteins\t$non_inteins\t$blastp_avg\n";
              }

              else{
                next;
              }
            }
            closedir BLASTP_EVAL;
          }
          closedir BLASTP;
        }

        #for psiblast results
        elsif($sample_contents eq "psiblast"){
          opendir PSIBLAST, "$directory\/$sample_number\/$sample_contents";
          while(my $psiblast_runs = readdir PSIBLAST){
            opendir PSIBLAST_EVAL, "$directory\/$sample_number\/$sample_contents\/$psiblast_runs";
            while(my $psiblast_sample = readdir PSIBLAST_EVAL){
              
              #directory of interest
              if($psiblast_sample eq "psiblast.blast"){
                #gets statistics
                my($psiblast_matches,$psiblast_inteins,$psiblast_avg)=PARSE_BLAST("$directory\/$sample_number\/$sample_contents\/$psiblast_runs\/$psiblast_sample");
                my $non_inteins = $psiblast_matches-$psiblast_inteins;

                #print to file
                print $psiblast_stats "$sample_number\t$psiblast_runs\t$psiblast_matches\t$psiblast_inteins\t$non_inteins\t$psiblast_avg\n";
              }

              else{
                next;
              }
            }
            closedir PSIBLAST_EVAL;
          }
          closedir PSIBLAST;
        
        
        }
        
        #skip otherwise
        else{
          next;
        }

      }
      closedir DIR_SUB;
    }

    else{
      next;
    }

  }
  closedir DIR;
}

sub PARSE_BLAST{
  #takes blast file as input
  #returns file information (#of matches, # of intein matches, average intein length)

  my $blast_file = shift;
  my %matches;
  my $num_matches = my $intein_matches = my $length_sum = 0;

  open(my $blast, "< $blast_file");
  while(<$blast>){
    chomp;
    my @match_split = split(/\t/,$_);
    #ignore secondary matches to the same database entry
    if($matches{$match_split[1]}){
      next;
    }
    else{
      $matches{$match_split[1]}="Y";
      $num_matches++;
      if($match_split[1]=~/intein/){
        $intein_matches++;
        $length_sum+= $match_split[3];
      }
      else{
        next;
      }
    }
  }
  close $blast;

  my $average_l = $length_sum/$intein_matches;

  return($num_matches,$intein_matches,$average_l);
}