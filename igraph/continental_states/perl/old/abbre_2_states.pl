#!/usr/bin/perl
# abbre_2_states.pl  

# First version: Aug 9, 2002 Hong Qin

# Preconditions:
#   Input file  1 : The old files containing the abbreviation for the states
#   Input file  2:  The translation table that for the states and their abbreviatons

# Postcoditions:
#   Output file 1 : The new file that containing the full state names

use strict; use warnings; 
use lib '/shar/lib/perl/'; # shanghai
use BeginPerlBioinfo;  use Util; use moran;
my $debug = 1;

if (! $ARGV[2]) { print "Usage perl $0 old_file translation_table new_file\n"; exit(0); }

my $infile1 = $ARGV[0];   my $infile2 = $ARGV[1];  my $outfile = $ARGV[2]; 

#get the old_file
my @old_lines = get_file_data( $infile1 );


# get the translation table
my %abbr_2_state = ();
my @lines = get_file_data( $infile2 );  
foreach my $line (@lines) {
  chomp $line;
   my ($state, $abbr, @rest) = split ( /[\s\t]+/, $line );
   $abbr = uc $abbr;
   $state = uc $state;
   if ( ( defined $state ) and ( defined $abbr ) ) {
     $abbr_2_state{$abbr} = $state;
   }
}

# translation starts
open ( OUT, ">$outfile");
foreach my $line (@old_lines) {
  chomp $line;
  my @elements = split ( /[\s\t\n]+/, $line);  
  my @new_elements = ();
  foreach my $element ( @elements) {
    my $tmp = uc $element;
    if ( exists $abbr_2_state{$tmp} ) {
	push ( @new_elements, $abbr_2_state{$tmp} );
    } else { 
	push ( @new_elements, $element );
    }
  }
  my $new_line = join ( "\t", @new_elements);
  print OUT "$new_line\n";
}

close (OUT);

exit;
