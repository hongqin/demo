#!/usr/bin/perl
# test_dijkasta.pl

# Nov 15, 2005 Hong Qin
# Oct 8, 2002 Hong Qin

use strict; use warnings; 
use Getopt::Long;
use lib '/home/hqin/lib/perl/'; # Shanghai Linux
use BeginPerlBioinfo;  use Util; use moran; use Graph;

my $network_file = undef;   
my $dic_file = undef;  
my $cat_file = undef;
my $color_file = undef;
my $null_report_file = undef;
my $out_network_file = undef; 
my $report_file = undef;
my $debug = 0;

if (! $ARGV[0]) { _help(); exit(1); }

GetOptions ('in_network=s' => \$network_file,   # required
            'out_network=s' => \$out_network_file,
            'report=s' => \$report_file,
	    'debug=i' => \$debug);
 
if ( ! $network_file) { _help(); exit(1); }
my $small_tmstamp =  get_short_time_stamp_US();   

# key data variables
 my %neighbours = ();		# the network in adj-list format
 my %dictionary = ();		# the network in adj-list format

# key variables for clustering
 
# key variables for statistics of the clustering result

# common temporary variables
 my ( @lines, $old_fasta_header, $fasta_header, $flag, %tmph ) =();
 my $time_stamp = undef;

# get the network configuration and the default dictionary
 $flag = parse_single_network_in_adjList(\$network_file, \%neighbours, \%dictionary, \$old_fasta_header);
 if ($debug>1) { print "\nDefault dictonary is:\n"; showHashTable(\%dictionary); } 

# here test the dijkastra subroutines
my $s = 'CALIFORNIA';
my %distances = ();
my %predecessors = ();
$flag = dijkstra(\%neighbours, $s, \%distances, \%predecessors );  
showHashTable(\%distances);
print "\n###################\n";
showHashTable(\%predecessors);

my $d = 'MAINE'; 
my $head = $d;
my $path = $d;
my $tail = '';
while ( $tail !~ $s ) {
  $tail = $predecessors{ $head };
  $path .= "<-$tail";
  $head = $tail;
}

print "\$path = $path\n";

open(TMP, ">__distances");
my @keys = sort(keys %distances);
hash_2_filehandle(*TMP, \%distances, \@keys, \"\t" );
close(TMP);


exit;

#
# end of main
#

#-----------------------------------------
sub _help {
  print "Usage:\n\tperl $0 \n\tRequired input fields:\n\t\t --in_network=flname1\n";
}

