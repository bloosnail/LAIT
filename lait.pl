#!/usr/bin/perl

###########################
#Daniel Hui
#dah124@pitt.edu
#University of Pittsburgh
#June 2015
###########################

use strict;
use warnings;

my $software = lc($ARGV[0]);

#LAMP
if($software eq "lamp"){
	if($ARGV[1] eq "2"){
        
        if(@ARGV != 5){
            print "USAGE:: perl lait.pl <lamp> <2> <map> <ped> <path>\n";
	    print "There are " . scalar @ARGV . " arguments instead of 5.\n";
            die;
            }	
       }elsif($ARGV[1] eq "3"){

	if(@ARGV != 8){
	    print "USAGE:: perl lait.pl <lamp> <3> <map> <ped> <freqs_pop1> <freqs_pop2> <freqs_pop3> <path>\n";
	    print "There are " . scalar @ARGV . " arguments instead of 8.\n";
		die;
		}
	}else{
        print "Sorry, LAIT does not currently support $ARGV[1]-way admixture for LAMP.\n";
        die;
        }

	print "Creating files...\n";

	open my $MAP, "$ARGV[2]", or die "mapfile open error";
        open my $PED, "$ARGV[3]", or die "pedfile open error";
	my $path = "$ARGV[-1]";
	
	my $directory = "$path";
  	unless(mkdir $directory) {
     		print "";
    	}
	print "Writing to $directory\n";
	
	open my $GENO, ">$path/geno.txt";
	open my $CONFIG, ">$path/config.txt";
	open my $SNPS, ">$path/chr.pos";

	my %snps;
	my @alts;
	
	if ($ARGV[1] eq "3"){
		open my $FREQ1, "$ARGV[4]" or die "freq1 error";
		open my $FREQ2, "$ARGV[5]" or die "freq2 error";
		open my $FREQ3, "$ARGV[6]" or die "freq3 error";
	
		open my $OUTFREQ1, ">$path/freqs_pop1.txt" or die "outfreq1 error";
		open my $OUTFREQ2, ">$path/freqs_pop2.txt" or die "outfreq2 error";
		open my $OUTFREQ3, ">$path/freqs_pop3.txt" or die "outfreq3 error";

		(%snps, @alts) = &lampfreq($FREQ1, $FREQ2, $FREQ3, $MAP, $OUTFREQ1, $OUTFREQ2, $OUTFREQ3, %snps, @alts);

#get alternate alleles from any frequency file (used pop1)
	seek($FREQ1, 0, 0);
	while(<$FREQ1>){
	  chomp;
	  my @line = split /\s+/;
	  if (exists $snps{$line[0]}){
	    push(@alts, $line[4]);
	    }
	  }
      }

	&lampconfig($CONFIG, $ARGV[1]);
	seek ($MAP, 0, 0) or die "can't go back: $!";

#lampsnps, didn't work when passed hash into sub
#"Odd number of elements in hash assignment"
if ($ARGV[1] eq "3"){

#SNP file
#lampsnps()
my $lc = 0;
my %counts;
while (<$MAP>){
        chomp;
        my @line = split(/\s+/, $_);
        if (exists $snps{$line[1]} && $line[3] ne "0"){
	$counts{$lc} = $lc;
        print $SNPS "$line[3]\n";
    }
	$lc++;
  }
	
	print "SNP file done.\n";

#Genotype file
#lampgeno()
while (<$PED>){
		my $count = 0;
                my @line = split /\s+/;
                #alt/variant is 0, reference is 1
                for (my $i = 0; $i < @line - 1; $i+=2){
	  	if (exists $counts{$i/2}){
                  if ( ($line[$i] eq $alts[$count]) && ($line[$i+1] eq $alts[$count]) ){
                        print $GENO "0 ";
                  }elsif ( ($line[$i] ne $alts[$count]) && ($line[$i+1] ne $alts[$count]) ) {
                        print $GENO "2 ";
                  }elsif ( ($line[$i] eq $alts[$count]) && ($line[$i+1] ne $alts[$count] )  ||  ($line[$i] ne $alts[$count]) && ($line[$i+1] eq $alts[$count] )  ) {
                        print $GENO "1 ";
                  }else{
                        print $GENO "-1 ";
                        }
			$count++;
                }
	     }
                print $GENO "\n";
        }
	print "Genotype file done.\n";

#2way
}else{
my %lines = &lampsnps($MAP, $SNPS);
&lampgeno($PED, $GENO, \%lines);
}

	print "All done! Check out $path\n";
	
        }



#LAMPLD
elsif($software eq "lampld" || $software eq "lamp-ld"){

my $path = "$ARGV[-1]";

	
	if($ARGV[1] eq "2"){
        
        if(@ARGV != 8){
            print "USAGE:: perl lait.pl <lampld> <2> <map> <ped> <snps> <haps_pop1> <haps_pop2> <path>\n";
	    print "There are " . scalar @ARGV . " arguments instead of 8.\n";
            die;
            }	

#3way
       } elsif($ARGV[1] eq "3"){

        print "Running LAMP-LD for 3way admixed samples.\n";

        if(@ARGV != 9){
            print "USAGE:: perl lait.pl <lamp-ld> <3> <map> <ped> <snps> <haps_pop1> <haps_pop2> <haps_pop3> <path>\n";
	    print "There are " . scalar @ARGV . " arguments instead of 9.\n";
            die;
            }
	    
  	    print "Creating files...\n";
        }else{
        print "Sorry, $ARGV[-1]-way admixture is not supported";
	die;
	}

        #END 3WAY
 
if($ARGV[1] eq "2"){
        print "Creating files...\n";
        }

	my $directory = "$path";
  	unless(mkdir $directory) {
     		print "";
    	}
	print "Writing to $directory\n";

	open my $MAP, "$ARGV[2]", or die "mapfile open error";
        open my $PED, "$ARGV[3]", or die "pedfile open error";
        open my $SNPS, "$ARGV[4]", or die "snps open error";
        open my $HAPS1, "$ARGV[5]", or die "pop1 haps open error";
        open my $HAPS2, "$ARGV[6]", or die "pop2 haps open error";
	
	open my $POS, ">$path/chr.pos";
	open my $GENO, ">$path/genofile.gen";
	open my $POP1, ">$path/pop1.hap", or die "pop1 geno open error";
        open my $POP2, ">$path/pop2.hap", or die "pop2 geno open error";
	
	my ($pedCount, $hapCount) = &ldpos($MAP, $SNPS, $POS);

	&ldgeno($PED, $GENO, $pedCount);
	&ldhap($HAPS1, $POP1, $hapCount);
	print "Reference haplotypes pop1 done.\n";
	&ldhap($HAPS2, $POP2, $hapCount);
	print "Reference haplotypes pop2 done.\n";
	
	if($ARGV[1] eq "3"){
            open my $HAPS3, "$ARGV[7]", or die "pop3 haps open error";
            open my $POP3, ">$path/pop3.hap", or die "pop3 geno error: $!";
            &ldhap($HAPS3, $POP3, $hapCount);
	    print "Reference haplotypes pop3 done.\n";
	}
	
	print "All done! Check out $path\n";
	
}

#ELAI
elsif($software eq "elai"){
	
	my $path = "$ARGV[-1]";
	
	#ELAI 2WAY
	if($ARGV[1] eq "2"){
	print "You have chosen ELAI for 2way admixed samples.\n";
	
	if(@ARGV != 9){
            print "USAGE:: perl lait.pl <elai> <2> <map> <ped> <haps_pop1> <snps_pop1> <haps_pop2> <snps_pop2> <path>\n";
	    print "There are " . scalar @ARGV . " arguments instead of 9.\n";
            die;
            }	
        }    
        #END 2WAY
        
        #ELAI 3WAY
        elsif($ARGV[1] eq "3"){
		
	print "You have chosen ELAI for 3way admixed samples.\n";	
		
	if(@ARGV != 11){
            print "USAGE:: perl lait.pl <elai> <3> <map> <ped> <haps_pop1> <snps_pop1> <haps_pop2> <snps_pop2> <haps_pop3> <snps_pop3> <path>\n";
	    print "There are " . scalar @ARGV . " arguments instead of 11.\n";
            die;
            }
        }
        #END 3WAY
        
	else{
            print "Sorry, LAIT does not currently support $ARGV[1]-way admixture.";		
            die;
	}
	
	#ELAI ANY WAY

	my $directory = "$path";
  	unless(mkdir $directory) {
     		print "";
    	}
	print "Writing to $directory\n";

	print "Creating files...\n";
	
	open my $MAP, "$ARGV[2]", or die "mapfile open error";
        open my $PED, "$ARGV[3]", or die "pedfile open error";
        open my $HAPS1, "$ARGV[4]", or die "pop1 haps open error";
        open my $SNPS1, "$ARGV[5]", or die "pop1 snps open error";
        open my $HAPS2, "$ARGV[6]", or die "pop2 haps open error";
        open my $SNPS2, "$ARGV[7]", or die "pop2 snps open error";
	
	open my $POS, ">$path/SNP.pos", or die "posfile open error";
        open my $GENO, ">$path/admix.geno", or die "genofile open error";
        open my $POP1, ">$path/pop1.geno", or die "pop1 geno open error";
        open my $POP2, ">$path/pop2.geno", or die "pop2 geno open error";
	
    &admixandpos($MAP, $PED, $POS, $GENO);
    print "Admixed genotype file done.";
    print "Position file done.";
    &sourcegeno($HAPS1, $SNPS1, $POP1);
    print "Reference genotype pop1 file done.";
    &sourcegeno($HAPS2, $SNPS2, $POP2);
    print "Reference genotype pop2 file done.";

       if($ARGV[1] eq "3"){
         open my $HAPS3, "$ARGV[8]", or die "pop3 haps open error";
         open my $SNPS3, "$ARGV[9]", or die "pop3 snps open error";
         open my $POP3, ">$path/pop3.geno", or die "pop3 geno open error";
         &sourcegeno($HAPS3, $SNPS3, $POP3);
	 print "Reference genotype pop3 file done.";
	}

    print "All done! Check out $path\n";

}


#HAPMIX
elsif($software eq "hapmix"){
	
		if($ARGV[1] eq "2"){
        
        if(@ARGV != 10){
            print "USAGE:: perl lait.pl <hapmix> <2> <chr#> <map> <ped> <ref_snps> <recombmap> <haps_pop1> <haps_pop2> <path>\n";
	    print "There are " . scalar @ARGV . " arguments instead of 10.\n";
            die;
            }	
        
        print "Creating files...\n";
        
	my $chrNo = $ARGV[2];
	open my $MAP, "$ARGV[3]", or die "mapfile open error";
        open my $PED, "$ARGV[4]", or die "pedfile open error";
        open my $SNPS, "$ARGV[5]", or die "snps open error";
        open my $RECMAP, "$ARGV[6]" or die "recmap open error";
        open my $HAPS1, "$ARGV[7]", or die "pop1 haps open error";
        open my $HAPS2, "$ARGV[8]", or die "pop2 haps open error";
	my $path = "$ARGV[-1]";
	
	my $directory = "$path";
  	unless(mkdir $directory) {
     		print "";
    	}
	print "Writing to $directory\n";

	my $dir2 = "$path/RUN";
  	unless(mkdir $dir2) {
     		print "$dir2 already exists\n";
    	}	
	
	open my $GENO, ">$path/AAgenofile.$chrNo" or die "genofile open error";
	open my $AASNPS, ">$path/AAsnpfile.$chrNo" or die "AAsnpfile open error";
	open my $POP1, ">$path/pop1genofile.$chrNo", or die "pop1 geno open error";
        open my $POP2, ">$path/pop2genofile.$chrNo", or die "pop2 geno open error";
        open my $SNPS1, ">$path/pop1snpfile.$chrNo", or die;
        open my $SNPS2, ">$path/pop2snpfile.$chrNo" or die;
        open my $INTER, ">$path/AAinterfile.txt" or die "interfile open error";
        open my $INTER1, ">$path/interfile1.txt" or die "interfile open error";
        open my $INTER2, ">$path/interfile2.txt" or die "interfile open error";
        open my $RATES, ">$path/rates.$chrNo" or die "rates open error";
	open my $PAR, ">$path/../parameters.$chrNo.par" or die ".par error";

	#get smallest and largest positions
	my $last;
	my $minPos;
	my $lc;
	while(<$RECMAP>){ 
	  $lc++;
	  $last = $_; 
	    if($lc == 2){
		my @fline = split /\s+/, $_;
		$minPos = $fline[0];	
	    }
	}
	my @line =  split /\s+/, $last;
	my $maxPos = $line[0];

	#make files
	&par($PAR, $chrNo, $path);
        &hapmixadmixsnps($MAP, $SNPS, $AASNPS, $SNPS1, $SNPS2, $maxPos, $minPos);
	print "AAsnpfile done.\n";	
	
	seek ($RECMAP, 0, 0) or die "can't go back";
	&hapmixrec($RECMAP, $RATES);
	print "Rates file done.\n";
	
	seek ($SNPS, 0, 0) or die "can't go back";
	seek ($MAP, 0, 0) or die "can't go back";
	
	
	&hapmixadmixgeno($SNPS, $MAP, $PED, $INTER, $GENO, $path, $maxPos, $minPos);
	print "AAgenofile done.\n";
	seek ($SNPS, 0, 0);
	seek ($MAP, 0, 0);
	&hapmixrefhaps($SNPS, $MAP, $HAPS1, $INTER1, $POP1, $path, "interfile1.txt", $maxPos, $minPos);
	print "Reference haps pop1 done.\n";

        seek ($SNPS, 0, 0);
        seek ($MAP, 0, 0);

	&hapmixrefhaps($SNPS, $MAP, $HAPS2, $INTER2, $POP2, $path, "interfile2.txt", $maxPos, $minPos);
	print "Reference haps pop2 done.\n";	

	print "All done! Check out $path\n";
	
	}else{
        print "HAPMIX can only do 2-way admixture.\n";       
	die;
	}

}

else{
    print "SEE INPUTS BY ENTERING:: perl lait.pl <software> <#waysAdmixed>\n";
    print "The software choices are lamp, lampld, elai, or hapmix.\n";
    print "The #waysAdmixed are 2 or 3.\n";
    die;
}




#<--------------------------#######SUBS########---------------------------->


#****************************LAMP****************************

################lamp freq
sub lampfreq{
	my ($FREQ1, $FREQ2, $FREQ3, $MAP, $OUT1, $OUT2, $OUT3, %all, @alts) = @_;

my %freq1;
my %freq2;
my %freq3;

while(<$FREQ1>){
	chomp;
	my @line = split /\s+/;
	$freq1{$line[0]}=$line[0];
}

while(<$FREQ2>){
	chomp;
	my @line = split /\s+/;
	$freq2{$line[0]}=$line[0];
}

while(<$FREQ3>){
	chomp;
	my @line = split /\s+/;
	$freq3{$line[0]}=$line[0];
}

#put each position from POSFILE into a hash as both the key and value
while (<$MAP>){
                chomp;
		my @line = split /\s+/;
		if( exists $freq1{$line[1]}) {
		if( exists $freq2{$line[1]}) {
		if( exists $freq3{$line[1]}) {
                      $all{$line[1]} = $line[1];
			}
		}
        }
}

#get the position value for each line in freqfile, and if it exists in positions hash then print out the 
#note that this is MINOR allele frequency

seek($FREQ1,0,0);
seek($FREQ2,0,0);
seek($FREQ3,0,0);

while (<$FREQ1>){
                chomp;
                my @line = split (/\s+/, $_);
                if (exists $all{$line[0]} )  {
                	print $OUT1 "$line[5]\n";
                }
        }
	print "Frequency file pop1 done.\n";

while (<$FREQ2>){
                chomp;
                my @line = split (/\s+/, $_);
                if (exists $all{$line[0]} )  {
                                print $OUT2 "$line[5]\n";
                }
        }
	print "Frequency file pop2 done.\n";

while (<$FREQ3>){
                chomp;
                my @line = split (/\s+/, $_);
                if (exists $all{$line[0]} )  {
                                print $OUT3 "$line[5]\n";
                }
        }
	print "Frequency file pop3 done.\n";

	return (%all);
}




################lamp snp list
sub lampsnps{
  my ($MAP, $SNPS) = @_;
  my $count = 0;
  my %lines;
  my %snps;
  while (<$MAP>){
    chomp;
    my @line = split /\s+/;

    if($line[3] ne "0"){
      if(exists $snps{$line[3]}){
        my $derp = 0;
        }else{
          $snps{$line[3]} = $line[3];
	  $lines{$count} = $count;  
          print $SNPS "$line[3]\n"; 
        }
      }
      $count++;
    }
print "SNP file done!\n";
return %lines;
}

################lamp ref geno
sub lampgeno{
my ($PED, $GENO, $lines_ref) = @_;
my %lines = %$lines_ref;

	while(<$PED>){
	  my @line = split(/\s+/, $_);
	  #print out 0, 1, or 2 depending on how many ref/alt alleles
		for (my $i = 0; $i < @line - 1; $i+=2){
	  	 if(exists $lines{$i/2}){
		  if ( ($line[$i] eq "A" || $line[$i] eq "T") && ($line[$i+1] eq "A"  || $line[$i+1] eq "T" ) ){
			print $GENO "2	";
		 }elsif ( ($line[$i] eq "C" || $line[$i] eq "G")  && ($line[$i+1] eq "C"  || $line[$i+1] eq "G" ) ) {
			print $GENO "0	";
		 }elsif ( ( ($line[$i] eq "A" || $line[$i] eq "T" ) && ($line[$i+1] eq "C" || $line[$i+1] eq "G") )   ||   ( ($line[$i] eq "C" || $line[$i] eq "G") && ($line[$i+1] eq "A" || $line[$i+1] eq "T" )  )  ) {
			print $GENO "1	";
		 }else{
			print $GENO "-1	";
			}	
		  }
	  	}
		print $GENO "\n";
             }
	print "Genotype file done.\n";
}


################lamp config file
sub lampconfig{
my ($CONFIG, $POPS) = @_;
					#no. of populations
print $CONFIG "# Number of populations\npopulations=$POPS\n\n"
."# To use LAMP with ancestral allele frequencies, provide files for allele frequencies for the pure populations\n";

#allele frequencies for the pure populations
if ($POPS == 3){ print $CONFIG "pfile=freqs_pop1.txt,freqs_pop2.txt,freqs_pop3.txt\n\n#######################################################################\n";}
else{ print $CONFIG "#pfile = freqs_pop1.txt,freqs_pop2.txt\n\n#######################################################################\n";}
print $CONFIG "# Required files\n#######################################################################\n"
		
		#genotype file			    #SNP pos file
."# Genotypes\ngenofile=geno.txt\n# SNP positions\nposfile=chr.pos\n# Output file of inferred ancestries.\n"

						#output file
."# Defaults to 'ancestry.txt'\noutputancestryfile=ancestry.txt\n\n"
."#######################################################################\n"
."# Parameters for the windowing scheme\n#######################################################################\n\n"

	#offset for adjacent windows				      #recombination rate
."# The offset of adjacent windows\noffset=0.2\n# Recombination rate\nrecombrate=1e-8\n";

	#no. of generations					#alpha
if ($POPS == 3) { print $CONFIG "# Number of generations\ngenerations=10\n# Alpha (must sum to one)\nalpha=0.7,0.1,0.2\n\n\n"}
else {
print $CONFIG "# Number of generations\ngenerations=7\n# Alpha (must sum to one)\nalpha=0.2,0.8\n\n\n"
}
print $CONFIG "#######################################################################\n"
."#######################################################################\n"

	#R^2 cutoff
."# R^2 Cutoff to  use for pruning SNPs\nldcutoff=0.1";
	print "Configuration file done.\n";	
}


#****************************LAMP-LD****************************


################lamp-ld ref haps
sub ldhap{
	
	my ($INHAP, $OUTHAP, $count_ref) = @_;
	my %count = %$count_ref;

	my @arr;
	for my $k (keys %count){
	push(@arr, $k);
}


	print scalar @arr."\n";
while (<$INHAP>){
		chomp;
		my @line = split("", $_);
		for (my $i = 0; $i < @line ; $i++){
		  if (exists $count{$i}){
			if ( $line[$i] eq "T" || $line[$i] eq "A"){ print $OUTHAP "0"; }
			elsif ( $line[$i] eq "G" || $line[$i] eq "C"){ print $OUTHAP "1"; } 
			else{ print $OUTHAP "?"; }
		  }
		}
		print $OUTHAP "\n";
	}

close $INHAP;
close $OUTHAP;
}


################lamp-ld snp pos
sub ldpos{
	my ($MAP, $SNPS, $POS) = @_;
	my %snps;
	my %subSnps;

	my %pc;
	my %hc;

	while(<$SNPS>){
		chomp;
		$snps{$_} = $_;
	}

#subset positions
#line counts for .ped file
	my $mc = 0;
	while (<$MAP>){
		chomp;
		my @line = split /\s+/;
	
	if( exists $snps{$line[1]} ){
	  if (exists $subSnps{$line[1]}){
	    my $derp;
	  }else{
		print $POS "$line[3]\n";
		$pc{$mc} = $mc;
		$subSnps{$line[1]} = $line[1];
		}			
	  $mc++;
	}
}

#line counts for hap file
	seek($SNPS, 0, 0) or die "can't go back: $!";

	my $sc=0;
	while (<$SNPS>){
	  chomp;
	  if(exists $subSnps{$_}){
	    $hc{$sc} = $sc;
	  }	
	    $sc++;
	}

	print "Position file done.\n";		
	return (\%pc, \%hc);
}


################lamp-ld ref geno
sub ldgeno{
	my ($PED, $GENO, $count_ref) = @_;
	
	my %count = %$count_ref;	

	while (<$PED>){
		chomp;
		my @line = split(/\s+/, $_);
		
		for (my $i = 0 ; $i < @line ; $i+=2){
		 if(exists $count{$i/2}){ 
		  if ("$line[$i]$line[$i+1]" eq "AA" || "$line[$i]$line[$i+1]" eq "AT" || "$line[$i]$line[$i+1]" eq "TA" || "$line[$i]$line[$i+1]" eq "TT"){
			print $GENO "0";
		  }elsif ("$line[$i]$line[$i+1]" eq "AG" || "$line[$i]$line[$i+1]" eq "GA" || "$line[$i]$line[$i+1]" eq "AC" || "$line[$i]$line[$i+1]" eq "CA" || "$line[$i]$line[$i+1]" eq "TG" || "$line[$i]$line[$i+1]" eq "GT" || "$line[$i]$line[$i+1]" eq "TC" || "$line[$i]$line[$i+1]" eq "CT"){
			print $GENO "1";
		  }elsif ("$line[$i]$line[$i+1]" eq "GG" || "$line[$i]$line[$i+1]" eq "GC" || "$line[$i]$line[$i+1]" eq "CG" || "$line[$i]$line[$i+1]" eq "CC"){
			print $GENO "2";
		  }else{
			print $GENO "?";
		  }
	   } 
	 }
	  print $GENO "\n"; 
      }
	print "Admixed genotype file done.\n";
    }

#****************************ELAI****************************

################elai, admix geno and snp pos
sub admixandpos{
	
my ($map, $ped, $pos, $geno) = @_;

my $individuals;
my $SNPs;
my @genoArray;

#lists no. of SNPs and rs no.
while (<$map>){
	chomp;
	my $line = $_;
	
	#pos
	my @mapArray = split(/\s+/, $line);
	print $pos "$mapArray[1] $mapArray[3] $mapArray[0]\n";
	
	#rs no.
	my @lineArray = split(/\s+/, $line);
	push (@genoArray, "$lineArray[1]");
	
	#no of SNPs
	$SNPs += 1;
}
	
#lists no. of individuals and geno
while (<$ped>){
	chomp;
	$individuals += 1;
	my @line = split(/\s+/, $_);

#rm headers if they're there
=pod
	for (my $i = 0; $i < 6; $i++){
		shift(@line);
	}
=cut
	#pair each 2 indices after that and delete old line
	my @newLine;
	my $lineCount=0;
	for (my $i = 0; $i < scalar @line; $i+=2){
			$newLine[$lineCount] = "$line[$i]$line[$i+1]";
			$lineCount++;
		}
	
	@line = ();
	
	#add new pairs to each SNP
	for (my $i = 0; $i < scalar @newLine; $i++){
		$genoArray[$i] = "$genoArray[$i] $newLine[$i]";
	}
}
	
#print statements
say $geno "$individuals ";
say $geno $SNPs;

foreach (@genoArray){
	print $geno "$_\n";
}

close $map;
close $ped;
close $geno;
}


################elai ref geno
sub sourcegeno{

my ($hap, $snp, $pop) = @_;

my $individuals;
my $SNPs;
my @genoArray;

#lists no. of SNPs and rs no.
while (<$snp>){
	chomp;
	push(@genoArray, $_);
	$SNPs += 1;
}
	
#lists no. of individuals and genotype
my $count = 0;
my @oldLine;

while (<$hap>){
	chomp;
	$individuals += 1;
	
	#do ops depending on every 2 lines
	if ($count == 0){
		@oldLine = split(//, $_);
	}
	
	if ($count == 1){
		my @line = split(//, $_);
	
		#add new pairs to each SNP
		for (my $i = 0; $i < @line; $i++){
			$genoArray[$i] = "$genoArray[$i],$oldLine[$i]$line[$i]";
		}
	}
	
	$count++;
	if ($count == 2){
			$count = 0;
		}
}

say $pop $individuals/2 ." =" ;
say $pop $SNPs;

foreach (@genoArray){
	print $pop "$_\n";
}

}

#****************************HAPMIX****************************
###############hapmix parameter
sub par{
	my ($PAR, $chrNo, $path) = @_;
	print $PAR
"GENOTYPE:1
OUTPUT_SITES:0
SITE_POSITIONS: 1 1000000000
THETA:0.2
LAMBDA:6.0
RECOMBINATION_VALS:600 900
MUTATION_VALS:0.2 0.2 0.01
MISCOPYING_VALS:0.05 0.05
REFPOP1GENOFILE:$path/pop1genofile.$chrNo
REFPOP2GENOFILE:$path/pop2genofile.$chrNo
REFPOP1SNPFILE:$path/pop1snpfile.$chrNo
REFPOP2SNPFILE:$path/pop2snpfile.$chrNo
ADMIXSNPFILE:$path/AAsnpfile.$chrNo
ADMIXGENOFILE:$path/AAgenofile.$chrNo
REF1LABEL:CEU
REF2LABEL:YRI
RATESFILE: $path/rates.$chrNo
ADMIXPOP: AA
CHR:$chrNo
OUTDIR:$path/RUN
HAPMIX_MODE:LOCAL_ANC
OUTPUT_DETAILS:ANC_INT_THRESH
THRESHOLD:0.0
KEEPINTFILES:0";

print "Parameter file done.\n";
}


################hapmix ref haps
sub hapmixrefhaps{
        my ($SNPS, $MAP, $HAP, $INTER, $OUTHAP, $path, $str, $maxPos, $minPos) = @_;

#read in and put admixed SNPs in hash
my %admixSNP;
while (<$MAP>){
                chomp;
		my @line = split /\t/;
	if ($line[3] < $maxPos && $line[3] > $minPos) {$admixSNP{$line[1]} = $line[1]};
        }

#cross admixed SNPs with those in simu .map 
my %crossSNPs;
my $pedCol = 0;

while(<$SNPS>){
                $pedCol++;
                chomp;

                if (exists($admixSNP{$_})){
                                $crossSNPs{$pedCol} = $pedCol;
                        }
        }

my $size = keys %crossSNPs;

#now that we have hash %crossSNPs w/ keys and values as no. of column of SNP that we want in .ped, 
#we skip each outline that is not the key/value that we want

while (<$HAP>) {

chomp;
my @line = split "";
#change letters to numbers
          for (my $i = 0; $i < @line ; $i++){

if (exists($crossSNPs{$i})){
                if ($line[$i] eq "A" || $line[$i] eq "T"){ print $INTER "0"; }
                elsif($line[$i] eq "G" || $line[$i] eq "C"){ print $INTER "1"; }
                else{ print $INTER "9"; }
        }

       }
print $INTER "\n";
}

close $INTER;
open $INTER, "<$path/$str" or die "error";

#check seek($INTER,0,0); instead^

my @outline;
my $lastcol;

while (<$INTER>){
        chomp;
  my @line = split("", $_);

      $lastcol = 0;
        $lastcol = $#line if $#line > $lastcol;
        my $oldlastcol = $lastcol;

        for (my $i = $oldlastcol; $i < $lastcol; $i++) {
                if (defined($oldlastcol)){
                $outline[$i] = "\t" x $oldlastcol;
                }
        }


        for (my $i=0; $i <=$lastcol; $i++) {
                $outline[$i] .= "$line[$i]\t"
        }
}

#rm whitespace and print
        for (my $i=0; $i <= $lastcol; $i++) {
                $outline[$i] =~ s/\s+//g;
                print $OUTHAP $outline[$i]."\n";
                                                            
}
close $INTER;
unlink "$path/$str" or print "can't delete ref interfile at $INTER: $!\n";

}


################hapmix admixed geno
sub hapmixadmixgeno{
	my ($SNPS, $MAP, $PED, $INTER, $GENO, $path, $maxPos, $minPos) = @_; 
	
#read in and put admixed SNPs in hash
my %admixSNP;
while (<$SNPS>){
                chomp;
                $admixSNP{$_} = $_;
        }

#cross admixed SNPs with those in simu .map 
my %crossSNPs;
my $pedCol = 0;

while(<$MAP>){
                $pedCol++;
                chomp;
                my @line = split /\t/;
              if ($line[3] < $maxPos && $line[3] > $minPos){
                if (exists($admixSNP{$line[1]})){
                                $crossSNPs{$pedCol} = $pedCol;
                        }
		}
        }

#now that we have hash %crossSNPs w/ keys and values as no. of column of SNP that we want in .ped, 
#we skip each outline that is not the key/value that we want

while (<$PED>) {
        
  chomp;
my @line = split /\s+/;

#change letters to numbers
        for (my $i = 0; $i < @line - 1; $i+=2){
                if (exists($crossSNPs{$i / 2})){
                #^ skips every 2 lines since pedigree file has two alleles per line
                         
                if ( ($line[$i] eq "A" && (($line[$i+1] eq  "A") || $line[$i+1] eq "T")) || ($line[$i] eq "T" && (($line[$i+1] eq  "A") || $line[$i+1] eq "T"))  ){
                        print $INTER "0";
                }

                elsif ( ($line[$i] eq "C" && (($line[$i+1] eq  "A") || $line[$i+1] eq "T")) || ($line[$i] eq "G" && (($line[$i+1] eq  "A") || $line[$i+1] eq "T"))
                || ($line[$i] eq "T" && (($line[$i+1] eq  "C") || $line[$i+1] eq "G")) || ($line[$i] eq "A" && (($line[$i+1] eq  "C") || $line[$i+1] eq "G")) ){
                        print $INTER "1";
                }

                elsif(  ($line[$i] eq "C" && (($line[$i+1] eq  "G") || $line[$i+1] eq "C")) || ($line[$i] eq "G" && (($line[$i+1] eq  "G") || $line[$i+1] eq "C") ) ){
                        print $INTER "2";
                }

                else{
                        print $INTER "9";
                }
        }
        }
        
        print $INTER "\n";
       }

close $INTER;
open $INTER, "$path/AAinterfile.txt" or die "erorr";

my @outline;
my $lastcol;

while (<$INTER>){
        chomp;
        my @line = split("", $_);
       
        $lastcol = 0;
        $lastcol = $#line if $#line > $lastcol;
        my $oldlastcol = $lastcol;

        for (my $i = $oldlastcol; $i < $lastcol; $i++) {
                $outline[$i] = "\t" x $oldlastcol;
        }


        for (my $i=0; $i <=$lastcol; $i++) {
                $outline[$i] .= "$line[$i]"
                }
        }


  for (my $i=0; $i <= $lastcol; $i++) {
                print $GENO $outline[$i]."\n";
}

close $INTER;
unlink "$path/AAinterfile.txt" or print "can't delete AA interfile\n";

}





################hapmix recombination rates
sub hapmixrec{
	
	my ($RECMAP, $RATES) = @_;
	
#read in RECMAPFILE, ignore first line
my $lineCount = 1;
my @outline;

while (<$RECMAP>) {
	
  if ($lineCount > 1){
  chomp;
  my @line = split(" ", $_);	

	#columns to rows
   my $lastcol = 0;
   $lastcol = $#line if $#line > $lastcol;
   my $oldlastcol = $lastcol;


        for (my $i = $oldlastcol; $i < $lastcol; $i++) {
                $outline[$i] = "\t" x $oldlastcol;
        }


        for (my $i=0; $i <=$lastcol; $i++) {
                $outline[$i] .= "$line[$i] " #change this char if want space or tab or whatever
                }
	
	}
	$lineCount++;
}

	print $RATES ":sites:";
	print $RATES ($lineCount - 2) . "\n";
        print $RATES $outline[0]."\n";
        print $RATES $outline[2]."\n";
}


################hapmix admixed snps
sub hapmixadmixsnps{
	my ($MAP, $SNPS, $OUT, $OUT1, $OUT2, $maxPos, $minPos) = @_;
	
#read in SNP file and put each value into %SNPhash, with key and value the same
my %SNPhash;
while (<$SNPS>){
		chomp;
		$SNPhash{$_} = $_;
	}

#test if SNP is in mapfile and if so, print
while (<$MAP>){
	chomp;
	my @lineArray = split(/\s+/, $_);
	if ($lineArray[3] < $maxPos && $lineArray[3] > $minPos){
		if (exists($SNPhash{$lineArray[1]})) {
	print $OUT "\t$lineArray[1]\t$lineArray[0]\t$lineArray[2]\t$lineArray[3]\n";
	print $OUT1 "\t$lineArray[1]\t$lineArray[0]\t$lineArray[2]\t$lineArray[3]\n";
	print $OUT2 "\t$lineArray[1]\t$lineArray[0]\t$lineArray[2]\t$lineArray[3]\n";
		}
		}
	}	
}
