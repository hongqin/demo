#!/usr/bin/perl
# adjList_2_adjMatrix.pl      

# First version: Aug 16, 2002 Hong Qin

# Preconditions:
#   Input file  1 : Network configuration file in the FASTA-liked, linked-list format

# Output:  adjacency matrix

use strict; use warnings; use lib '/shar/lib/perl/'; use BeginPerlBioinfo;  use Util; use moran;
my $debug = 2;

if (! $ARGV[1]) { print "Usage perl $0 adjacency_list.config output\n"; exit(0); }
my $infile1 = $ARGV[0];		my $outfile = $ARGV[1]; 

my ( @lines, $line );  # input buffer
my ( $i, $j, $k ) = ();
my ( %category_to_index, @categories ) = ();  #category lookup tables

#get the network configuration
@lines = get_file_data( $infile1 );
my $old_fasta_header = shift @lines;  # this is a non-intelligent approach!

# set up %neighbours
 my %neighbours = ( );
# setup %redundant_pairs (including both id1-id2 and id2-id1 pairs )
my %redundant_pairs = ();
foreach my $line ( @lines) {
  my @current_neighbours = ();
  chomp $line;
  my ( $head, @rest ) = split ( /[\s\t]+/, $line );
  my ( $head_id, $head_category ) = split ( /:/, $head );
  foreach my $partner ( @rest ) {
    my ( $id, $category ) = split ( /:/, $partner );
    if ( defined $id ) {  push ( @current_neighbours, $id );   }
    $neighbours{$head_id} = join ( ' ', @current_neighbours );
    $redundant_pairs{ $head.'  '.$id } ++;
    $redundant_pairs{ $id.'  '.$head } ++;
  }
}
if ($debug==2) { print " Network is:\n"; showHashTable(\%neighbours); }

# setup @ordered_ids
my @ordered_ids =  sort ( keys %neighbours );

# output the 2d array
open (OUT, ">$outfile");
print OUT "@ordered_ids \n";
foreach my $i ( 0 .. $#ordered_ids ) {
  foreach my $j  ( 0 .. $#ordered_ids ) {
    my $edge = $ordered_ids[$i] .'  '. $ordered_ids[$j] ;
    if (exists $redundant_pairs{$edge}) { print OUT "1  " }
    else { print OUT "0  "; }
  }
  print OUT "\n";
}
close (OUT);

exit;
