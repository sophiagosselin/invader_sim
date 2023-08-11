#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;
use Cwd qw(cwd);

#this will get me some statistics

my @evals=("1e-3","1e-5","1e-8","1e-10","1e-13","1e-15","1e-18","1e-20");
my @ids=("0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95");

MAIN();


sub MAIN{
    my $home_dir = cwd;

    foreach my $evalue (@evals){
        my $evaldir = $evalue;
        $evaldir=~s/\-/\_/g;
        
        foreach my $identity (@ids){
            #create subdirectories
            my($iddir)=($identity=~/0\.(.*)/);

            mkdir("ice_blast_runs/$iddir\_$evaldir");

            #copy files
            opendir ICE, "ice_blast_runs";
            while (my $files = readdir ICE){
                
                #skip over directories
                if(-d "$files"){
                    next;
                }

                #work only on the files
                else{
                    copy("ice_blast_runs\/$files","ice_blast_runs/$iddir\_$evaldir\/$files");
                }

            }
            closedir ICE;

            chdir "ice_blast_runs/$iddir\_$evaldir";
            system("perl iceblast.pl -in random_intein.fasta -psidb intein_subset.fasta -outdb all_sample_sequences.fasta -t 10 -id $identity -e $evalue -ds .25");
            chdir "$home_dir";
        }
    }   
}