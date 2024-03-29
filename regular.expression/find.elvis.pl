#! /usr/bin/perl -w

BEGIN { unshift(@INC,"/home/hqin/lib/perl/");   }

use strict; use warnings; 
use Getopt::Long;   
use Util;  use BeginPerlBioinfo; use FASTA;

 my $in_fasta_fl = "prot.db.faa"; 
 my $out_fasta_fl = "_finding.elvis.faa";	

 my $tmp_fasta_fl = "/tmp/_tmp.demo.fasta"; 
 my $flag = paddle_fasta_file1($in_fasta_fl, $tmp_fasta_fl );

 open (IN, "<$tmp_fasta_fl" ); my @lines = <IN>; close (IN); 

 my $count = 0;
 open (OUT, ">$out_fasta_fl");
 for( my $i=1; $i<=$#lines; $i++ ) {
    #if ($lines[$i] =~ /ELVISISFRMMEMPHIS/) {
    #if ($lines[$i] =~ /ELVISIS.*MEMPHIS/) { # Pattern search using regular expression
    if ($lines[$i] =~ /ELVIS[IS|CAME].*MEMPHIS/) { # Pattern search using regular expression
      	print OUT $lines[$i-1]; 	#print out FASTA header
	print OUT $lines[$i];  		#print out sequence
	$count ++; 			#tally the counts
    }
 }
  close (IN);
  close (OUT);

 print "\t Pattern is found $count times. \n"; 

 system (  "rm $tmp_fasta_fl " );

exit;

