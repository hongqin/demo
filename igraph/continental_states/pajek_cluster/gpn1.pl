#!/usr/bin/perl
# generate_pajek_network_b0.1.pl 

# Upated on Aug 12, 02 Hong Qin
#	Now parse the standard network file in adjacency-list

# Updated on July 22, 02  Hong Qin
#    New feature to read network in adjacency-list and clustering output from color_clusters.pl

# Version b0.0 : July 11, 2002 Hong Qin 
#	Only parse the pairwise data

# Preconditions:
#   Use *.ctl file as input which indicates all the options and data files
#  required fields:
#	OUTPUT_FILE = some_file_name.net
#  	NETWORK_FORMAT = pairwise OR adjacency_list
#	NETWORK_FILE  = some_file_name
#  optional fields:
#	CLUSTER_FLAG = YES or NO
#	DICTIONARY_FILE  = NONE or some_file_name
#	REMOVE_EDGES_NOT_IN_DICTIONARY = YES or NO
#	CATEGORY_TO_COLOR_FILE = NONE or some_file_name
#	COLOR_EDGE_BY_PAIRING_NODES = YES or NO
#	CONSIDER_SELF_INTERACTION = YES or NO
#	..., ... etc.

use strict; use warnings; 
# use lib '/shar/lib/perl/';   	# Linux
 use lib '\lib\perl/';  	# Windows
use Pajek;	use Util;	use moran;
my $debug = 1;

# temporary variables
my ( $line, $tmp, @lines, @input_lines, $i, @elements, $field, $input ) =();
my ( $id1, $id2, @rest) = ();
my ( $arc, $pair, $edge ) = ();
my ( $category, $color)= (); 

# key variables
my $ctl_file = "gp.ctl";  	if ( $ARGV[0]) {  $ctl_file = $ARGV[0]  }  # the control file
my $default_vertex_size = 1;
my ( %config ) = ();		# the configuration hash
my ( %nr_ids, %nr_pairs ) = ();		# the non_redundant nodes and edges

my ( @pajek_vertices, %positions, %pajek_edges, %pajek_arcs, @colors, @categories ) = ();
# @pajek_vertices is the array of vertices
# %positions stores the indice for the vertices in the array
# %pajek_edges are for undirected edges
# %pojek_arcs are for directed edges

my ( %dictionary, %input_dictionary ) = ();  # the dictionaries
my ( %cluster_sizes ) = ();		# the cluster_size hash
my ( %color_by_category, %category_by_color ) = ();    # the coloring scheme

open (DEBUG, ">_gpn_debug.txt" );

unless ( open ( CTL, "<$ctl_file") ) { die "Cannot open the configuration file ($ctl_file).\n" ; }  
@lines = <CTL>;  close (CTL);
@input_lines = grep ( /=/, @lines );
foreach $line ( @input_lines ) {
    ($field, $input) = split ( /=/, $line );
    $field =~ s/[\s\t\n]+//g;  	$input =~ s/[\s\t\n]+//g;  $field = uc $field;
    $config{$field} = $input;
    if ($debug) { print DEBUG "Field:[$field]\tValue:[$input]\n"; }
}


open ( IN, "<$config{'NETWORK_FILE'}"  ); 
@lines = <IN>; close (IN);
if ($debug) { print DEBUG "\n There are $#lines records in $config{'NETWORK_FILE'} "; }
if ($config{'NETWORK_FORMAT'} eq 'pairwise' ) {
 foreach $line (@lines) {
   ($id1, $id2, @rest ) = split ( /\s+/, $line );
   if ( $config{"CONSIDER_SELF_INTERACTION"} eq 'NO' ) {
     if ( ( $id1 gt '' ) and ( $id2 gt '') and ( $id1 ne $id2 ) ) {
        my ( $big, $small ) = order_big2small ( $id1, $id2 );
        $nr_ids{ $big } ++;  $nr_ids{ $small } ++;
        $nr_pairs{ $big."\t".$small } ++;
     }
   } else { # $config{"CONSIDER_SELF_INTERACTION"} eq 'YES' 
     if ( ( $id1 gt '' ) and ( $id2 gt '') ) {
        my ( $big, $small ) = order_big2small ( $id1, $id2 );
        $nr_ids{ $big } ++;  $nr_ids{ $small } ++;
        $nr_pairs{ $big."\t".$small } ++;       
     }
   } 
 }#$line loop
}#pairwise if loop
elsif (($config{'NETWORK_FORMAT'} eq 'adjacency_list') and ($config{'CLUSTER_FLAG'} eq 'YES') ) {
 shift @lines; # skip the fasta header

 foreach $line (@lines) {
  my ( @current_colors, @current_cluster_heads ) = ();
  my @clusters = split (  /\t+/, $line );

  foreach my $cluster (@clusters) {
    my ($cluster_list, $color) = split( /:/, $cluster);
    $cluster_list =~ s/\(//g;
    $cluster_list =~ s/\)//g;
    $color =~ s/[\s\t\n]+//g;
    @elements = split( /\s+/, $cluster_list);
    push ( @current_colors, $color );
    if ( defined $elements[1] ) {
      push ( @current_cluster_heads, $elements[0]);
      $nr_ids{$elements[0]} ++;
      $cluster_sizes{$elements[0]} = scalar @elements;
      $input_dictionary{$elements[0]} = $color;
    } else { 
      push ( @current_cluster_heads, $cluster_list) ; # single node cluster 
      $nr_ids{$cluster_list} ++;
      $cluster_sizes{$cluster_list} = 1;
      $input_dictionary{$cluster_list} = $color;
    }
  } # $cluster

  # now generate edges among the current clusters
  my $head = shift @current_cluster_heads;
  my $head_color = shift @current_colors;
  foreach my $node (@current_cluster_heads) {
    my ($big, $small) = order_big2small($head, $node);
    $nr_pairs{$big."\t".$small} ++;
  }

 } # $line loop
} # cluster if loop
elsif (($config{'NETWORK_FORMAT'} eq 'adjacency_list') and ($config{'CLUSTER_FLAG'} =~ /NO/i) ) {
 shift @lines; # skip the fasta header

 foreach $line (@lines) {
  my ( $head, @rest ) = split ( /\t/, $line );  # node are tab delimited
  my ( $head_id, $head_category ) = split ( /:/, $head );
  $head_category =~ s/\s+//g;
  $head_id =~ s/\s+//g;
  if ( $head_id gt '') {  $nr_ids{$head_id} ++;  }
  if (($head_id gt '') and ($head_category gt '')) {
     $input_dictionary{$head_id} = $head_category;
  }
  foreach my $partner ( @rest ) {
    my ( $id, $category ) = split ( /:/, $partner );
    my ($big, $small) = order_big2small($head_id, $id);
    $nr_pairs{$big."\t".$small} ++;
  } 
 } # $line loop
} # adjacency_list loop if ($debug) { print "\nThe input dictionary:\n"; showHashTable(\%input_dictionary) }

@pajek_vertices = keys %nr_ids;		# the vertices are generated here

foreach $i ( 0 .. $#pajek_vertices ) {
   $positions{$pajek_vertices[$i]} = $i;
}

foreach $pair ( keys %nr_pairs ) {
  ( $id1, $id2, @rest ) = split ( /\t/, $pair );
  # $arc = ( $positions{$id1} +1) ."  ". ( $positions{$id2} +1 );
  # $pajek_arcs{ $arc } ++;
  $edge = ( $positions{$id1} +1) ."  ". ( $positions{$id2} +1 );
  $pajek_edges{$edge} ++;
}

if ( ( $config{'DICTIONARY_FILE'} ne 'NONE' ) and (  $config{'DICTIONARY_FILE'} ne '' ) ) {
  open ( IN, "<$config{'DICTIONARY_FILE'}" ); @lines = <IN>; close (IN);
  foreach $line ( @lines ) {
    ( $id1, $category, @rest ) = split ( /[\s\t]+/, $line );
    if ( $id1 ne '' ) {
       if ( $category ne '' ) {
          $dictionary{$id1} = $category;
          $color_by_category{$category} = 'default';  # ??
       } else {
          $dictionary{$id1} = 'none';
       }
    } 
  }  
} else {
 foreach my $key ( keys %input_dictionary ) {
   $dictionary{$key} = $input_dictionary{$key};
 }
}


if ( ( $config{'CATEGORY_TO_COLOR_FILE'} ne 'NONE' ) and (  $config{'CATEGORY_TO_COLOR_FILE'} ne '' ) ) { 
  open ( IN, "<$config{'CATEGORY_TO_COLOR_FILE'}" ); @lines = <IN>; close (IN);
  foreach $line ( @lines ) {
    ( $category, $color, @rest ) = split ( /[\s\t]+/, $line );
    if ( $category ne '' ) {
       $color_by_category{$category} = $color;
       $category_by_color{$color} = $category;
    }
  }
} else {
  # Automatically generate a coloring scheme of the vertices
  @colors =  get_defined_Pajek_colors();
  @categories = keys %color_by_category;
  foreach $i ( 0 .. $#categories ) {
    $i = $i % 78;
    $color_by_category{ $categories[$i] } = $colors[$i];
  }
}

if ($debug) { print "Category_2_color translation table:\n"; showHashTable(\%color_by_category);}

# output
open ( OUT, ">$config{'OUTPUT_FILE'}" ); 

# output @pajek_vertices first
print OUT "*Vertices  " . ( $#pajek_vertices + 1) ." \n";
foreach $i ( 0 .. $#pajek_vertices ) {
  if ( $config{'CLUSTER_FLAG'} eq 'YES' ) {
     $tmp = '"' . $pajek_vertices[$i]. " ($cluster_sizes{$pajek_vertices[$i]})" . '"';
#     my $vertex_size = $default_vertex_size + round2(log($cluster_sizes{$pajek_vertices[$i]})/log(4) , 1);
	 my $vertex_size = round2(  sqrt( $cluster_sizes{$pajek_vertices[$i]} ) , 1);
     print OUT   "  ". ($i+1). "  $tmp ";
     print OUT " x_fact $vertex_size y_fact $vertex_size ";
  } else {
     $tmp = '"' . $pajek_vertices[$i] . '"'  ;   # how to put quote here ?????
     print OUT   "  ". ($i+1). "  $tmp ";
  }
  if ( (keys %dictionary ) > 2 )  {
     $category =  $dictionary{ $pajek_vertices[$i] };
     if ( (defined $category) and ( exists $color_by_category{ $category  } ) ) {
          print OUT " ic $color_by_category{ $category  } ";
     } else {
          print OUT " ic Gray ";    # default is Gray !!!
     }
  }
  print OUT " \n";
}

# output edges
print OUT "*Edges \n";
my ($n1, $n2, $category1, $category2, $color1, $color2, $edge_color) = ();
foreach $edge ( keys %pajek_edges ) {
  if ( ($config{'DICTIONARY_FILE'} ne 'NONE') and (exists $config{'COLOR_EDGE_BY_PAIRING_NODES'})
        and ( $config{'COLOR_EDGE_BY_PAIRING_NODES'} eq 'YES' )) {
     ( $n1, $n2 ) = split ( /[\s\t]+/, $arc );
     $n1 --; $n2 --;
     $category1 = $dictionary{ $pajek_vertices[$n1] };
     $category2 = $dictionary{ $pajek_vertices[$n2] };
     if( (defined $category1) and (defined $category2)
         and (exists $color_by_category{$category1} ) and (exists $color_by_category{$category2} )) {
         $color1 = $color_by_category{$category1};
         $color2 = $color_by_category{$category2};
         if ( $color1 eq $color2 ) { $edge_color = $color1 }
         else { $edge_color = 'Magenta' }  # this should be changed later
     }  
  } else { # the default color for edges
    $edge_color = 'Black';			# default edge color
  }
  print OUT "  $edge 1 c $edge_color \n";
}

close (OUT);
close (DEBUG);
exit;
