#!/usr/bin/perl -w
use strict;
use warnings;

my $directory = "ice_blast_runs";

opendir DIR, $directory ;
while(my $sample_num = readdir DIR){
  # $dir will be the sample number
  if(-d "$directory\/$sample_num"){
    opendir SUBDIR, "$directory\/$sample_num";
    while(my $contents = readdir SUBDIR){
      if($contents eq "blastp"){
        system("rm -r $directory\/$sample_num\/$contents\n");
      }
      elsif($contents eq "psiblast"){
        system("rm -r $directory\/$sample_num\/$contents\n");
      }
      else{
        next;
      }
    }
  }
}
