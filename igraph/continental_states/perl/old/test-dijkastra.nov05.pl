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
            'category_file=s' => \$cat_file,	
            'dictionary|d=s' => \$dic_file,
            'color_scheme=s' => \$color_file,
            'out_network=s' => \$out_network_file,
            'null_report=s' => \$null_report_file,
            'report=s' => \$report_file,
	    'debug=i' => \$debug);
 
if ( ! $network_file) { _help(); exit(1); }
my $small_tmstamp =  get_short_time_stamp_US();   
#if ( ! $out_network_file ) {  $out_network_file = '_cclu_'.$network_file.".".$small_tmstamp; }
#if ( ! $report_file) { $report_file = '_ccluReport_'.$network_file.".".$small_tmstamp; } 

if ($debug) { # show the fields
  print "Inputs:\n\t\$network_file= $network_file\n";
  if ($cat_file ) { print "\t\$cat_file= '$cat_file'\n";}
  else { print "\t\$cat_file= undef\n";}
  if (defined $dic_file) { print "\t\$dic_file= '$dic_file'\n"; }
  else { print "\t\$dic_file = undef\n"; }
  if (defined $color_file) { print "\t\$color_file= '$color_file'\n"; }
  else { print "\t\$color_file = undef\n"; }
  if (defined $null_report_file) { print "\t\$null_report_file= '$null_report_file'\n"; }
  else { print "\t\$null_report_file = undef\n"; }
  print "\t\$out_network_file= '$out_network_file'\n";
  print "\t\$report_file= '$report_file'\n";
  print "\t\$debug= $debug\n";
}

# key data variables
 my %neighbours = ();		# the network in adj-list format
 my %dictionary = ();		# the default dictionary in $network_file
 my %in_dictionary = ();	# the optional dictionary to override the default one
 my %category_to_color = (); 	# the optional color-scheme
 my @categories = ();		# the ordered categories for output
 my %category_2_index = ();	# hash each category to a ordered position
 my @in_categories = (); 	# the optional ordered categories

# key variables for clustering
 my %same_color_clusters = ( );  # store clustering result for each node, 
 my %nr_clusters = ();           # the unique clusters for final output 
 my %node_2_cluster_head = ();   # lookup table for node-to-cluster conversion
 my %neighboured_clusters = ();  # the final network in clusters, adjacency_list format
 my %new_pairs = ();  # This the new pairwise interaction hash, values are occurence for each pair
 		      # The pair is ordered from big2small
 
# key variables for statistics of the clustering result
 my $cluster_num = undef;  # total number of colored_clusters
 my %cluster_sizes = ();	 	# store the size for each cluster
 my %cluster_sizes_by_color = ();	# store the sizes for each color (category)
 my %largest_clusterSize_by_color = ();
 my %average_clusterSize_by_color = ();
 my %clusterNum_by_color = ();

# key varialbes for significance analysis
 # see below

# key variables for parsing and analyzing null models
 my %null_maxSizes_by_color = ();  	# string_container
 my %null_avgSizes_by_color = ();	# string_container
 my %null_clusterNums_by_color = ();	# string_container
 my $num_of_nulls = undef;		# number of null models

# common temporary variables
 my ( @lines, $old_fasta_header, $fasta_header, $flag, %tmph ) =();
 my $time_stamp = undef;

# get the optional @in_categories
 if (defined $cat_file) {
   $flag = parse_ordered_ids_by_row(\$cat_file,\@in_categories);  
 }

# get the optional coloring_scheme file
 if ( defined $color_file ) {
  $flag = parse_dictionary(\$color_file, \%category_to_color);
  if ($debug) { print "Color_scheme is:\n"; showHashTable(\%category_to_color); }
 }

# get the optional dictionary
 if (defined $dic_file) {
   if (defined $color_file ) {
     $flag = parse_dictionary(\$dic_file, \%in_dictionary, \%category_to_color);
   } else {
      $flag = parse_dictionary(\$dic_file, \%in_dictionary);
   }
   if ($debug>1) { print "\nOptional dictonary is:\n"; showHashTable(\%in_dictionary); } 
 }


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

