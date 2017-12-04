#!/usr/bin/perl
# classify_states.pl   

# First version: Aug 9, 2002 Hong Qin

# Preconditions:
#   Input file  1 : File contain the states and the official starting year

# Postcoditions:
#   Output file 1 : File that label each state with a category (initial, intermediate, recent)

use strict; use warnings; 
use lib '/shar/lib/perl/'; # shanghai
use BeginPerlBioinfo;  use Util; use moran;
my $debug = 0;

if (! $ARGV[1]) { print "Usage perl $0 input output \n"; exit(0); }

my $infile = $ARGV[0];  my $outfile = $ARGV[1]; 

#get the in_file
my @in_lines = get_file_data( $infile );

# classification starts
open ( OUT, ">$outfile");
foreach my $line (@in_lines) {
    chomp $line;
    my $category = '';
    my ( $state, $year, @rest ) = split ( /[\s\t\n]+/, $line);  
    if ( $year <= 1790 ) { $category = 'initial' }
    elsif (( $year >= 1781) and ( $year <= 1850 )) { $category = 'intermediate' }
    else  { $category = 'recent' }
    if ($debug) { print "[$state]\t[$year]\t[$category]\n";}
    print OUT "$state\t$category\n";
}
close (OUT);

exit;
