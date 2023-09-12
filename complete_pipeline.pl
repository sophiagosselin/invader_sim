#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use Cwd;


#globals
my @evals=("1e-3","1e-5","1e-8","1e-10","1e-13","1e-15","1e-18","1e-20");
my @ids=("0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95");
my $home_dir = cwd;
my $threads = $ARGV[0];

MAIN();

sub MAIN{

    #run invader sim w/o subsampling - edit desired parameters here!
    INVADER_SIM("-sp 1000 -b 1 -d 1 -sf .1 -m .1\n-sp 100 -b 2 -d 1 -sf .5 -m 1\n-sn 100 -ws 100 -evo evolver\n-nn 500 -a .5 -cat 4 -m 1\n-nn 300 -a 1 -cat 4 -m 1");

    #----ICEBLAST baseline experiment----#
    #copy all files into each sample's directory
    my($hashref1,$hashref2,$hashref3)=BASELINE_EXPERIMENT("baseline_original_experiment");

    #----Unified samples experiment----#
    #make unified sample fasta file
    my($hashref4,$hashref5,$hashref6)=UNIFIED_EXPERIMENT("baseline_original_experiment","baseline_unified_experiment");

    #----Cleanup----#
    mkdir("invader_sim_no_sub_sampling");
    rmove("simulated_sequences","invader_sim_no_sub_sampling\/simulated_sequences");
    rmove("invaded_sequences","invader_sim_no_sub_sampling\/invaded_sequences");
    rmove("extein","invader_sim_no_sub_sampling\/extein");
    rmove("intein","invader_sim_no_sub_sampling\/intein");

    #----Subsample Experiment----#
    #run invader sim w/ subsampling - edit desired parameters here!
    INVADER_SIM("-sp 1000 -b 1 -d 1 -sf .1 -m .1\n-sp 1000 -b 2 -d 1 -sf .5 -m 1\n-sn 100 -ws 100 -evo evolver -sub 100\n-nn 500 -a .5 -cat 4 -m 1\n-nn 300 -a 1 -cat 4 -m 1");
    
    #----Baseline---#
    my($hashref7,$hashref8,$hashref9)=BASELINE_EXPERIMENT("subsample_original_experiment");
    
    #----Unified----#
    my($hashref10,$hashref11,$hashref12)=UNIFIED_EXPERIMENT("subsample_original_experiment","subsample_unified_experiment");

    #----Cleanup 2----#
    mkdir("invader_sim_with_sub_sampling");
    rmove("simulated_sequences","invader_sim_with_sub_sampling\/simulated_sequences");
    rmove("invaded_sequences","invader_sim_with_sub_sampling\/invaded_sequences");
    rmove("extein","invader_sim_with_sub_sampling\/extein");
    rmove("intein","invader_sim_with_sub_sampling\/intein");

    #----Statistics----#
    STATISTICS($hashref1,$hashref2,$hashref3,$hashref4,$hashref5,$hashref6,$hashref7,$hashref8,$hashref9,$hashref10,$hashref11,$hashref12);

}

sub STATISTICS{
    my %original_baseline_blastp_files = %{ my $hr1 = shift};
    my %original_baseline_psiblast_files = %{ my $hr2 = shift};
    my %original_baseline_iceblast_files = %{ my $hr3 = shift};
    my %original_unified_blastp_files = %{ my $hr4 = shift};
    my %original_unified_psiblast_files = %{ my $hr5 = shift};
    my %original_unified_iceblast_files = %{ my $hr6 = shift};
    my %subsample_baseline_blastp_files = %{ my $hr7 = shift};
    my %subsample_baseline_psiblast_files = %{ my $hr8 = shift};
    my %subsample_baseline_iceblast_files = %{ my $hr9 = shift};
    my %subsample_unified_blastp_files = %{ my $hr10 = shift};
    my %subsample_unified_psiblast_files = %{ my $hr11 = shift};
    my %subsample_unified_iceblast_files = %{ my $hr12 = shift};

    #original sample blastp
    PARSE_BLAST_HASH(\%original_baseline_blastp_files,"original_baseline_blastp_stats.txt",0);

    #original sample psiblast
    PARSE_BLAST_HASH(\%original_baseline_psiblast_files,"original_baseline_psiblast_stats.txt",0);

    #original sample iceblast
    PARSE_FASTA_HASH(\%original_baseline_iceblast_files,"original_baseline_iceblast_stats.txt",0);

    #unified sample blastp
    PARSE_BLAST_HASH(\%original_unified_blastp_files,"unified_baseline_blastp_stats.txt",1);

    #unified sample psiblast
    PARSE_BLAST_HASH(\%original_unified_psiblast_files,"unified_baseline_psiblast_stats.txt",1);

    #unified sample iceblast
    PARSE_FASTA_HASH(\%original_unified_iceblast_files,"unified_baseline_iceblast_stats.txt",1);

    #subsample blastp
    PARSE_BLAST_HASH(\%subsample_baseline_blastp_files,"subsample_baseline_blastp_stats.txt",0);

    #subsample psiblast
    PARSE_BLAST_HASH(\%subsample_baseline_psiblast_files,"subsample_baseline_psiblast_stats.txt",0);

    #subsample iceblast
    PARSE_FASTA_HASH(\%subsample_baseline_iceblast_files,"subsample_baseline_iceblast_stats.txt",0);

    #unified subsample blastp
    PARSE_BLAST_HASH(\%subsample_unified_blastp_files,"subsample_unified_blastp_stats.txt",1);

    #unified subsample psiblast
    PARSE_BLAST_HASH(\%subsample_unified_psiblast_files,"subsample_unified_psiblast_stats.txt",1);

    ##unified subsample iceblast
    PARSE_FASTA_HASH(\%subsample_unified_iceblast_files,"subsample_unified_iceblast_stats.txt",1);

}

sub PARSE_FASTA_HASH{

    my %hash_of_fasta_results = %{my $ref = shift};
    my $output_file = shift;
    my $unified_toggle = shift;

    open(my $out, "+> $output_file");
    print $out "Sample#\tEvalue\tIdentity\tTotal_Matches\tIntein_Matches\tFalse_Positives\tAverage_Intein_Length\n";
    
    if($unified_toggle == 1){
        foreach my $evalue (keys %hash_of_fasta_results){
            foreach my $identity (keys %{$hash_of_fasta_results{"$evalue"}}){
                my $resultsfile = $hash_of_fasta_results{"$evalue"}{"$identity"};
                my($total_matches,$int_matches,$avg_intein_l)=PARSE_FASTA_OUTPUT($resultsfile);
                my $false_positives = ($total_matches-$int_matches);
                print $out "1\t$evalue\t$identity\t$total_matches\t$int_matches\t$false_positives\t$avg_intein_l\n";
            }
        }
    }

    else{
        foreach my $sample (keys %hash_of_fasta_results){
            foreach my $evalue (keys %{$hash_of_fasta_results{"$sample"}}){
                foreach my $identity (keys %{$hash_of_fasta_results{"$sample"}{"$evalue"}}){
                    my $resultsfile = $hash_of_fasta_results{"$sample"}{"$evalue"}{"$identity"};
                    my($total_matches,$int_matches,$avg_intein_l)=PARSE_FASTA_OUTPUT($resultsfile);
                    my $false_positives = ($total_matches-$int_matches);
                    print $out "$sample\t$evalue\t$identity\t$total_matches\t$int_matches\t$false_positives\t$avg_intein_l\n";
                }
            }
        }
    }
    close $out;
}

sub PARSE_FASTA_OUTPUT{
    #takes fasta as input
    #returns number of intein matches, total matches, and average intein length
    my $results_fasta = shift;
    my $num_matches = my $intein_matches = my $length_sum = 0;

    my %output_sequences=READIN_FASTA("$results_fasta");
    
    foreach my $asc (keys %output_sequences){
        $num_matches++;
        if($asc=~/w_intein/){
            $intein_matches++;
            $length_sum+=length($output_sequences{$asc});
        }
        else{
            next;
        }
    }

    my $average_l;
    if($num_matches == 0){
        $average_l = 0;
    }
    else{
        $average_l = $length_sum/$num_matches;
    }
      
    return($num_matches,$intein_matches,$average_l);
}

sub PARSE_BLAST_HASH{   
    my %hash_of_blast_results = %{my $ref = shift};
    my $output_file = shift;
    my $unified_toggle = shift;

    open(my $out, "+> $output_file");
    print $out "Sample#\tEvalue\tTotal_Matches\tIntein_Matches\tFalse_Positives\tAverage_Intein_Length\n";

    #check if toggle is triggered
    if($unified_toggle == 1){
        foreach my $evalue (keys %hash_of_blast_results){
            my $resultsfile = $hash_of_blast_results{"$evalue"};
            my($total_matches,$int_matches,$avg_intein_l)=PARSE_BLAST_OUTPUT($resultsfile);
            my $false_positives = ($total_matches-$int_matches);
            print $out "1\t$evalue\t$total_matches\t$int_matches\t$avg_intein_l\n";
        }
    }

    else{
        foreach my $sample (keys %hash_of_blast_results){
            foreach my $evalue (keys %{$hash_of_blast_results{"$sample"}}){
                my $resultsfile = $hash_of_blast_results{"$sample"}{"$evalue"};
                my($total_matches,$int_matches,$avg_intein_l)=PARSE_BLAST_OUTPUT($resultsfile);
                my $false_positives = ($total_matches-$int_matches);
                print $out "$sample\t$evalue\t$total_matches\t$int_matches\t$avg_intein_l\n";
            }
        }
    }
    close $out;
}
    


sub PARSE_BLAST_OUTPUT{
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

    my $average_l;
    if($intein_matches == 0){
        $average_l = 0;
    }
    else{
        $average_l = $length_sum/$intein_matches;
    }

    return($num_matches,$intein_matches,$average_l);
}

sub BLASTP{
    #note that output is in blastp.blast
    #returns file location
    my $query = shift;
    my $outdb = shift;
    my $evalue = shift;
    system("blastp -db $outdb -query $query -evalue $evalue -out blastp.blast -outfmt 6");

    my $path = cwd;
    return("$path\/blastp.blast");
}

sub PSIBLAST{
    #note that the output from the run will be found in psiblast.blast
    my $query = shift;
    my $psidb = shift;
    my $outdb = shift;
    my $evalue = shift;
    
    system("psiblast -db $psidb -query $query -out temp.txt -out_pssm intein.pssm -inclusion_ethresh $evalue -outfmt \"6 sseqid sstart send\" -num_iterations 4 -save_pssm_after_last_round");
    system("psiblast -db $outdb -in_pssm intein.pssm -out psiblast.blast -inclusion_ethresh $evalue -evalue $evalue -outfmt 6");

    my $path = cwd;
    return("$path\/psiblast.blast");
}

sub UNIFIED_EXPERIMENT{
    #runs ice-blast experiment on a unified dataset
    #additionally runs psiblast and blastp on the unified dataset
    my $origin_directory = shift;
    my $exp_directory = shift;
    my $counter = -1;
    my(%blastp_r,%psiblast_r,%iceblast_r);

    mkdir("$exp_directory");
    my $unique_counter = 0;
    open(my $unified, "+> $exp_directory\/all_sample_sequences.fasta");
    my @invaded_seq_files = glob "invaded_sequences/*";
    foreach my $invaded_fasta (@invaded_seq_files){
        my(%sequence_data) = READIN_FASTA($invaded_fasta);
        foreach my $asc (keys %sequence_data){
            print $unified "$asc\_$unique_counter\n";
            print $unified "$sequence_data{$asc}\n";
            $unique_counter++;
        }
        $counter++;
    }
    close $unified;

    #pick random intein to use as query and copy files over
    my $random_query = int(rand($counter))+1;

    copy("$origin_directory\/$random_query\/random_intein.fasta","$exp_directory\/random_intein.fasta");
    copy("$origin_directory\/$random_query\/intein_subset.fasta","$exp_directory\/intein_subset.fasta");
    
    #record this random selection
    open(my $log, "+> $exp_directory\/log.txt");
    print $log "Random intein is from sample $random_query\n";
    close $log;
    
    #makeblastdb
    MAKEBLASTDB("$exp_directory\/intein_subset.fasta");
    MAKEBLASTDB("$exp_directory\/all_sample_sequences.fasta");
    my @unified_files = glob "$exp_directory\/*";

    #directory setup
    mkdir("$exp_directory\/psiblast");
    mkdir("$exp_directory\/blastp");

    #paramater files
    foreach my $evalue (@evals){
        my $evaldir = $evalue;
        $evaldir=~s/\-/\_/g;
        
        #blastp / psiblast setup
        mkdir("$exp_directory\/psiblast\/$evaldir");
        mkdir("$exp_directory\/blastp\/$evaldir");

        #mkdirs copy files
        foreach my $identity (@ids){
            my($iddir)=($identity=~/0\.(.*)/);
            mkdir("$exp_directory\/$iddir\_$evaldir"); 
            foreach my $file (@unified_files){
                my($filename)=($file=~/$exp_directory\/(.*)/);
                copy($file,"$exp_directory\/$iddir\_$evaldir\/$filename");
                copy($file,"$exp_directory\/psiblast\/$evaldir\/$filename");
                copy($file,"$exp_directory\/blastp\/$evaldir\/$filename");
            }

            #run ICEBLAST for given experiment.
            copy("iceblast.pl","$exp_directory\/$iddir\_$evaldir\/iceblast.pl");
            chdir "$exp_directory\/$iddir\_$evaldir";
            my($result)=ICEBLAST("random_intein.fasta","intein_subset.fasta","all_sample_sequences.fasta",$identity,$evalue);
            chdir "$home_dir";
            $iceblast_r{"$evalue"}{"$identity"}=$result;
        }

        #blastp + psiblast exps.
        chdir "$exp_directory\/psiblast\/$evaldir";
        my($psiresult)= PSIBLAST("random_intein.fasta","intein_subset.fasta","all_sample_sequences.fasta",$evalue);
        chdir "$home_dir";
        $psiblast_r{"$evalue"}=$psiresult;
        
        chdir "$exp_directory\/blastp\/$evaldir";
        my($blastpresult)= BLASTP("random_intein.fasta","all_sample_sequences.fasta",$evalue);
        chdir "$home_dir";
        $blastp_r{"$evalue"}=$blastpresult;
    }

    return(\%blastp_r,\%psiblast_r,\%iceblast_r);
}

sub BASELINE_EXPERIMENT{
    #runs ice-blast experiment with no modifications (see the unified experiment sub)
    #additionally runs the psiblast and blastp experiments on the baseline files.
    my $exp_directory = shift;
    my(%blastp_r,%psiblast_r,%iceblast_r);

    mkdir("$exp_directory");
    my $counter = 1;
    my %all_sample_data;
    my @invaded_seq_files = glob "invaded_sequences/*";
    foreach my $invaded_fasta (@invaded_seq_files){
        mkdir("$exp_directory\/$counter");
        my($hashref_sampledata)=COPY_FASTAS($invaded_fasta,"$exp_directory\/$counter",$counter);
        %all_sample_data=MERGE_NESTED_HASHES($hashref_sampledata,\%all_sample_data);
        $counter++;
    }

    #create crucial files w/in each directory for baseline experiment
    opendir my $samples_home_dir, "$exp_directory";
    while(my $indv_sample_dir = readdir $samples_home_dir){
        if(-d "$exp_directory\/$indv_sample_dir"){
            next if("$indv_sample_dir" eq "." || "$indv_sample_dir" eq "..");
            
            #create databases and queries
            CREATE_INTEIN_DERIVATIVES(\%all_sample_data,$indv_sample_dir,"$exp_directory\/$indv_sample_dir");
            MAKEBLASTDB("$exp_directory\/$indv_sample_dir\/intein_subset.fasta");
            MAKEBLASTDB("$exp_directory\/$indv_sample_dir\/invaded_seqs.fasta");
            my @sample_files = glob "$exp_directory\/$indv_sample_dir\/*";

            #directory setup
            mkdir("$exp_directory\/$indv_sample_dir\/psiblast");
            mkdir("$exp_directory\/$indv_sample_dir\/blastp");

            #create + populate experimental parameter folders
            foreach my $evalue (@evals){
                my $evaldir = $evalue;
                $evaldir=~s/\-/\_/g;

                #blastp / psiblast setup
                mkdir("$exp_directory\/$indv_sample_dir\/psiblast\/$evaldir");
                mkdir("$exp_directory\/$indv_sample_dir\/blastp\/$evaldir");

                foreach my $identity (@ids){
                    my($iddir)=($identity=~/0\.(.*)/);
                    mkdir("$exp_directory\/$indv_sample_dir\/$iddir\_$evaldir"); 
                    foreach my $file (@sample_files){
                        my($filename)=($file=~/baseline\_experiment\/$indv_sample_dir\/(.*)/);
                        copy($file,"$exp_directory\/$indv_sample_dir\/$iddir\_$evaldir\/$filename");
                        copy($file,"$exp_directory\/$indv_sample_dir\/psiblast\/$evaldir\/$filename");
                        copy($file,"$exp_directory\/$indv_sample_dir\/blastp\/$evaldir\/$filename");
                    }

                    #run ICEBLAST for given experiment.
                    copy("iceblast.pl","$exp_directory\/$indv_sample_dir\/$iddir\_$evaldir\/iceblast.pl");
                    chdir "$exp_directory\/$indv_sample_dir\/$iddir\_$evaldir";
                    my($result)=ICEBLAST("random_intein.fasta","intein_subset.fasta","invaded_seqs.fasta",$identity,$evalue);
                    chdir "$home_dir";
                    $iceblast_r{"$indv_sample_dir"}{"$evalue"}{"$identity"}=$result;
                }

                #blastp + psiblast exps.
                chdir "$exp_directory\/$indv_sample_dir\/psiblast\/$evaldir";
                my($psiresult)= PSIBLAST("random_intein.fasta","intein_subset.fasta","invaded_seqs.fasta",$evalue);
                chdir "$home_dir";
                $psiblast_r{"$indv_sample_dir"}{"$evalue"}=$psiresult;
        
                chdir "$exp_directory\/$indv_sample_dir\/blastp\/$evaldir";
                my($blastpresult)= BLASTP("random_intein.fasta","invaded_seqs.fasta",$evalue);
                chdir "$home_dir";
                $blastp_r{"$indv_sample_dir"}{"$evalue"}=$blastpresult;

            }
        }
        else{
            next;
        }
    }
    closedir $samples_home_dir;

    return(\%blastp_r,\%psiblast_r,\%iceblast_r);
}

sub ICEBLAST{
    #returns path to output file

    my $query = shift;
    my $psidb = shift;
    my $outdb = shift;
    my $id = shift;
    my $ev = shift;
    system("perl iceblast.pl -in $query -psidb $psidb -outdb $outdb -t $threads -id $id -e $ev -ds .25");

    my $path = cwd;
    return("$path\/output\/all_matches.fasta");
}

sub MAKEBLASTDB{
    #does what it says on the tin

    my $input_file = shift;
    system("makeblastdb -in $input_file -dbtype prot -parse_seqids -out $input_file");
}

sub MERGE_NESTED_HASHES{
    #read the tin

	my @hashes_to_merge = @_;
	my %new_hash;
	foreach my $hashref (@hashes_to_merge){
		my %temp_hash = %{$hashref};
		%new_hash = (%new_hash,%temp_hash);
	}
	my $test = (\%new_hash);
	return(%new_hash);
}

sub CREATE_INTEIN_DERIVATIVES{
    #creates file containing only 1 random intein
    #creates file containing subset of inteins sans the random one
    #returns name of randomly selected intein

    my %file_data = %{my $hashref = shift};
    my $sample_number = shift;
    my $directory = shift;

    my %intein_data = READIN_FASTA($file_data{$sample_number}{"intein"});

    #create fasta file w/ randomly selected intein as query
    my @ascs = keys %intein_data;
    my $random_key = splice(@ascs,rand(@ascs),1);
    open(my $rand_fasta, "+> $directory\/random_intein.fasta"); 
    print $rand_fasta "$random_key\n";
    print $rand_fasta "$intein_data{$random_key}\n";
    close $rand_fasta;

    #create fasta file w/ subset of randomly selected inteins sans the query one above
    my $num_ascs = @ascs;
    my $prop_of_ascs = $num_ascs*0.1;
    open(my $subset_fasta, "+> $directory\/intein_subset.fasta");
    for(my $rand_counter=0; $rand_counter<=$prop_of_ascs; $rand_counter++){
        my $random_asc = splice(@ascs,rand(@ascs),1);
        print $subset_fasta "$random_asc\n";
        print $subset_fasta "$intein_data{$random_asc}\n";
    }
    close $subset_fasta;

}

sub INVADER_SIM{
    #creates param file from string then runs runs invader sim
    my $parameter_string = shift;

    open(my $param, "+> is.param");
    print $param $parameter_string;
    close $param;

    system("perl invader_sim.pl");
}

sub READIN_FASTA{
    #also standardizes any file it reads in
    my $infile = shift;
    my $accession="";
    my %sequences;
    
    open(IN, "< $infile");
    while(<IN>){
        chomp;
        $_=~s/[\/\-]/\_/g;
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

sub COPY_FASTAS{
    #takes an invaded sequence file output from invader sim
    #and a directory to copy into
    #finds all files associated with that experiment
    #copies files into directory.
    #returns nested hash of what files are in the directory

    my $target_experiment = shift;
    my $directory = shift;
    my $sample_num = shift;
    my %nested_hash;
    my($extein_num,$intein_num)=($target_experiment=~/invaded_sequences\/extein\_sample\_(.*?)\_intein\_sample\_(.*?)\.fasta/);

    my $extein_fasta = "simulated_sequences\/extein\/$extein_num\/sample_$extein_num.fasta";
    my $intein_fasta = "simulated_sequences\/intein\/$intein_num\/sample_$intein_num.fasta";

    copy($extein_fasta,"$directory\/extein_sample_$extein_num.fasta");
    copy($intein_fasta,"$directory\/intein_sample_$intein_num.fasta");
    copy($target_experiment,"$directory\/invaded_seqs.fasta");
    
    $nested_hash{$sample_num}{"extein"}="$directory\/extein_sample_$extein_num.fasta";
    $nested_hash{$sample_num}{"intein"}="$directory\/intein_sample_$intein_num.fasta";
    $nested_hash{$sample_num}{"invaded"}="$directory\/invaded_seqs.fasta";

    return(\%nested_hash);
}