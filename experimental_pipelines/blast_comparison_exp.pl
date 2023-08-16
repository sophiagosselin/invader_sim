#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;
use Cwd qw(cwd);

my $directory = "ice_blast_runs";
my @evals=("1e-3","1e-5","1e-8","1e-10","1e-13","1e-15","1e-18","1e-20");

MAIN();


sub MAIN{
  opendir DIR, $directory ;
  while(my $sample_num = readdir DIR){
    # $dir will be the sample number
    if(-d "$directory\/$sample_num"){
      next if($sample_num eq "." || $sample_num eq "..");
      mkdir("$directory\/$sample_num\/psiblast");
      mkdir("$directory\/$sample_num\/blastp");

      my $home_dir = cwd;
      foreach my $evalue (@evals){
        #create subdirectories
        my $dir_name = $evalue;
        $dir_name=~s/\-/\_/g;
        mkdir("$directory\/$sample_num\/psiblast\/$dir_name");
        mkdir("$directory\/$sample_num\/blastp\/$dir_name");

        #copy files
        opendir DIR_SUB, "$directory\/$sample_num";
        while (my $sample_contents = readdir DIR_SUB){
          #skip over directories
          if(-d "$directory\/$sample_num\/$sample_contents"){
            next;
          }
          #work only on the files
          else{
            copy("$directory\/$sample_num\/$sample_contents","$directory\/$sample_num\/psiblast\/$dir_name");
            copy("$directory\/$sample_num\/$sample_contents","$directory\/$sample_num\/blastp\/$dir_name");
          }
        }

        #run BLASTP
        chdir "$home_dir\/$directory\/$sample_num\/blastp\/$dir_name";
        #print "Test $home_dir\/$directory\/$sample_num\/blastp\/$dir_name\n";
        #system("export BLASTDB=$home_dir\/$directory\/$sample_num\/blastp\/$dir_name");
        system("blastp -db invaded_seqs.fasta -query random_intein.fasta -evalue $evalue -out blastp.blast -outfmt 6");
        chdir "$home_dir";

        #run PSIBLAST
        chdir "$home_dir\/$directory\/$sample_num\/psiblast\/$dir_name";
        #system("export BLASTDB=$home_dir\/$directory\/$sample_num\/psiblast\/$dir_name");
        system("psiblast -db intein_subset.fasta -query random_intein.fasta -out temp.txt -out_pssm intein.pssm -inclusion_ethresh $evalue -outfmt \"6 sseqid sstart send\" -num_iterations 4 -save_pssm_after_last_round");
        system("psiblast -db invaded_seqs.fasta -in_pssm intein.pssm -out psiblast.blast -inclusion_ethresh $evalue -evalue $evalue -outfmt 6");
        chdir "$home_dir";
      }
    }
    else{
      next;
    }
  }
}
