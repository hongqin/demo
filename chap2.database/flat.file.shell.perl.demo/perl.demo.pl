#!/usr/bin/perl
# demo
# pick out the genes expressed in liver

use strict; use warnings;

my $fl = "db.txt";

my %genes = ();

open (IN, "<$fl");
while (my $line =<IN>) {
 if ($line =~ /liver/ ) {
   my ($id, $tissue, @rest) = split( /\t|\s+/, $line);
   $genes{ $id } = $tissue;
 }
}
close (IN);

my $count = 1;
foreach my $gene (sort (keys %genes) ) {
 print "  $count  $gene\n";
 $count ++;
}
exit;


