#!/usr/bin/perl
# test_dijkasta.pl

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
if ( ! $out_network_file ) {  $out_network_file = '_cclu_'.$network_file.".".$small_tmstamp; }
if ( ! $report_file) { $report_file = '_ccluReport_'.$network_file.".".$small_tmstamp; } 

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
open(TMP, ">__distances");
my @keys = sort(keys %distances);
hash_2_filehandle(*TMP, \%distances, \@keys, \"\t" );
close(TMP);

my %dists_by_cats = ();
$flag = single_source_shortestDistances_by_categories(\%neighbours, $s, \%dictionary, \%dists_by_cats );
showHashTable(\%dists_by_cats);

my %avg_dists_by_cats = ();
$flag = all_sources_shortestDistances_by_categories(\%neighbours, \%dictionary, \%avg_dists_by_cats );    
showHashTable(\%avg_dists_by_cats);

# override the default dictionary if indicated so 
 if (defined $dic_file) { 
   %dictionary = %in_dictionary ; 
   if ($debug) { print "\nDefault dictionary is overridden.\n"; }
 }

# set the default @categories, %category_2_index by %dictionary
 foreach my $definiton (values %dictionary) {
   if (defined $definiton ) {
     $tmph{$definiton} ++;
   } else {
     print "\n\tNull definiton found. Bye.\n"; 
     exit(1);
   }
 }
 @categories = sort( keys %tmph );  # get the default oder for output

# override the default @categories by the optional input category-order
 if (defined $cat_file) {
   foreach my $cat (@in_categories) {
     if ( !(exists $tmph{$cat}) ) {  # consistency check
       print "\n\tCategories in '$cat_file' are incompatible. Bye.\n";
       exit(1);
     }
   }
   @categories = @in_categories;  # override
   if ($debug) { print "\nUse the input category order for output.\n"; }
 }
 if ($debug) {print "Ordered categories for output are:[@categories]\n"}

 foreach my $i ( 0..$#categories ) {
   $category_2_index{$categories[$i]} = $i;
 }


#-------here do the same-color clustering-----------

# find the clusters of the same_colored_neighbours for each node
 $flag = all_sources_BFS_on_same_colored_nodes(\%neighbours,\%same_color_clusters,\%dictionary ); 
 if ($debug==2) { showHashTable(\%same_color_clusters); }

 open ( DEBUG, ">_debug_color_cluster");

# Decompose the graph into clusters

 # map each to the largest symbol in its cluster
 foreach my $node (keys %same_color_clusters ) {
   my @elements = split ( /[\s\t\n]+/, $same_color_clusters{$node}  );
   my $largest_symbol =  get_largest_symbol(@elements); 
   if ( ! exists $nr_clusters{$largest_symbol} ) {  
     $nr_clusters{ $largest_symbol } = $same_color_clusters{$node}  ;
   }
   $node_2_cluster_head{$node} =  $largest_symbol;   
 }

 # Replace the all nodes in %neighbours with sets of cluster_heads
 $flag = translate_adjacency_list_2_pairs_by_lookupTable (\%neighbours, \%new_pairs, \%node_2_cluster_head );
 if ($debug) { print "\nThe redundant pairwise interactions:\n"; showHashTable(\%new_pairs); }

 %new_pairs = remove_homodimers(\%new_pairs);
 if ($debug) { print "\nThe unique pairwise interactions:\n"; showHashTable(\%new_pairs); }







exit;

#
# end of main
#

#-----------------------------------------
sub _help {
  print "Usage:\n\tperl $0 \n\tRequired input fields:\n\t\t --in_network=flname1\n";
  print "\tOptional input fields:\n";
  print "\t\t --dictionary=flname3 \n";
  print "\t\t --category_file=flname2\n";
  print "\t\t --color_scheme=flname3 \n";
  print "\t\t --null_report= output from 'ccluall' \n";
  print "\t\t --debug=value (0, 1, 2, default is 0)\n";
  print "\tOptional output fields:\n";
  print "\t\t --out_network=flname3 (default: _cclu_.'\$network_file'.)\n";
  print "\t\t --report=some_flname (default: _ccluReport_.'\$network_file'.)\n";
  print "\t\t --tmp_file=flname4 (default: _debug_cclu_.'\$network_file'.)\n";
  return;
}

#-----------------------------------------
# Usage: cclu_set_hash_from_line(\%hash, $line);
sub cclu_set_hash_from_line {
 use strict; use warnings; 
 my ($ref_hash, $line) = @_;

 my (@elements) = split ( /,/,  $line);
 foreach my $element (@elements) {
     $element =~ s/\s+//g;
     my($cat,$value,@rest) = split ( /=>/, $element);
     if ( (defined $cat) and (defined $value) ) {
      $$ref_hash{$cat} .= $value . ' ';
     }
 }
 return;
}

