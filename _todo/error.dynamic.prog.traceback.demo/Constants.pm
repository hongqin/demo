# The DNA nucleotides, the amino acids, and the standard genetic code.
# Here encapsulated in a package.
# KVL Biolinux * sestoft@itu.dk * 2002

use strict;
use warnings;
 
@Constants::nucs = split(//, "ACGT");
@Constants::aas  = split(//, "ARNDCQEGHILKMFPSTWYV");
           
%Constants::gencode = 
    (TTT => "F", TTC => "F", TTA => "L", TTG => "L",
     TCT => "S", TCC => "S", TCA => "S", TCG => "S",
     TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
     TGT => "C", TGC => "C", TGA => "*", TGG => "W",
     
     CTT => "L", CTC => "L", CTA => "L", CTG => "L",
     CCT => "P", CCC => "P", CCA => "P", CCG => "P", 
     CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
     CGT => "R", CGC => "R", CGA => "R", CGG => "R", 
     
     ATT => "I", ATC => "I", ATA => "I", ATG => "M",
     ACT => "T", ACC => "T", ACA => "T", ACG => "T", 
     AAT => "N", AAC => "N", AAA => "K", AAG => "K",
     AGT => "S", AGC => "S", AGA => "R", AGG => "R",

     GTT => "V", GTC => "V", GTA => "V", GTG => "V", 
     GCT => "A", GCC => "A", GCA => "A", GCG => "A",
     GAT => "D", GAC => "D", GAA => "E", GAG => "E",
     GGT => "G", GGC => "G", GGA => "G", GGG => "G" );

# The BLOSUM45 amino acid substitution matrix

@Constants::blosum45 = 
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
  ( [  5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-2,-2, 0],  # A
    [ -2, 7, 0,-1,-3, 1, 0,-2, 0,-3,-2, 3,-1,-2,-2,-1,-1,-2,-1,-2],  # R
    [ -1, 0, 6, 2,-2, 0, 0, 0, 1,-2,-3, 0,-2,-2,-2, 1, 0,-4,-2,-3],  # N
    [ -2,-1, 2, 7,-3, 0, 2,-1, 0,-4,-3, 0,-3,-4,-1, 0,-1,-4,-2,-3],  # D
    [ -1,-3,-2,-3,12,-3,-3,-3,-3,-3,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1],  # C
    [ -1, 1, 0, 0,-3, 6, 2,-2, 1,-2,-2, 1, 0,-4,-1, 0,-1,-2,-1,-3],  # Q
    [ -1, 0, 0, 2,-3, 2, 6,-2, 0,-3,-2, 1,-2,-3, 0, 0,-1,-3,-2,-3],  # E
    [  0,-2, 0,-1,-3,-2,-2, 7,-2,-4,-3,-2,-2,-3,-2, 0,-2,-2,-3,-3],  # G
    [ -2, 0, 1, 0,-3, 1, 0,-2,10,-3,-2,-1, 0,-2,-2,-1,-2,-3, 2,-3],  # H
    [ -1,-3,-2,-4,-3,-2,-3,-4,-3, 5, 2,-3, 2, 0,-2,-2,-1,-2, 0, 3],  # I
    [ -1,-2,-3,-3,-2,-2,-2,-3,-2, 2, 5,-3, 2, 1,-3,-3,-1,-2, 0, 1],  # L
    [ -1, 3, 0, 0,-3, 1, 1,-2,-1,-3,-3, 5,-1,-3,-1,-1,-1,-2,-1,-2],  # K
    [ -1,-1,-2,-3,-2, 0,-2,-2, 0, 2, 2,-1, 6, 0,-2,-2,-1,-2, 0, 1],  # M
    [ -2,-2,-2,-4,-2,-4,-3,-3,-2, 0, 1,-3, 0, 8,-3,-2,-1, 1, 3, 0],  # F
    [ -1,-2,-2,-1,-4,-1, 0,-2,-2,-2,-3,-1,-2,-3, 9,-1,-1,-3,-3,-3],  # P
    [  1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-3,-1,-2,-2,-1, 4, 2,-4,-2,-1],  # S
    [  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1, 2, 5,-3,-1, 0],  # T
    [ -2,-2,-4,-4,-5,-2,-3,-2,-3,-2,-2,-2,-2, 1,-3,-4,-3,15, 3,-3],  # W
    [ -2,-1,-2,-2,-3,-1,-2,-3, 2, 0, 0,-1, 0, 3,-3,-2,-1, 3, 8,-1],  # Y
    [  0,-2,-3,-3,-1,-3,-3,-3,-3, 3, 1,-2, 1, 0,-3,-1, 0,-3,-1, 5]   # V 
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
    );

# The BLOSUM50 amino acid substitution matrix

@Constants::blosum50 = 
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
  ( [  5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-2,-1,-1,-3,-1, 1, 0,-3,-2, 0],  # A
    [ -2, 7,-1,-2,-4, 1, 0,-3, 0,-4,-3, 3,-2,-3,-3,-1,-1,-3,-1,-3],  # R
    [ -1,-1, 7, 2,-2, 0, 0, 0, 1,-3,-4, 0,-2,-4,-2, 1, 0,-4,-2,-3],  # N
    [ -2,-2, 2, 8,-4, 0, 2,-1,-1,-4,-4,-1,-4,-5,-1, 0,-1,-5,-3,-4],  # D
    [ -1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1],  # C
    [ -1, 1, 0, 0,-3, 7, 2,-2, 1,-3,-2, 2, 0,-4,-1, 0,-1,-1,-1,-3],  # Q
    [ -1, 0, 0, 2,-3, 2, 6,-3, 0,-4,-3, 1,-2,-3,-1,-1,-1,-3,-2,-3],  # E
    [  0,-3, 0,-1,-3,-2,-3, 8,-2,-4,-4,-2,-3,-4,-2, 0,-2,-3,-3,-4],  # G
    [ -2, 0, 1,-1,-3, 1, 0,-2,10,-4,-3, 0,-1,-1,-2,-1,-2,-3, 2,-4],  # H
    [ -1,-4,-3,-4,-2,-3,-4,-4,-4, 5, 2,-3, 2, 0,-3,-3,-1,-3,-1, 4],  # I
    [ -2,-3,-4,-4,-2,-2,-3,-4,-3, 2, 5,-3, 3, 1,-4,-3,-1,-2,-1, 1],  # L
    [ -1, 3, 0,-1,-3, 2, 1,-2, 0,-3,-3, 6,-2,-4,-1, 0,-1,-3,-2,-3],  # K
    [ -1,-2,-2,-4,-2, 0,-2,-3,-1, 2, 3,-2, 7, 0,-3,-2,-1,-1, 0, 1],  # M
    [ -3,-3,-4,-5,-2,-4,-3,-4,-1, 0, 1,-4, 0, 8,-4,-3,-2, 1, 4,-1],  # F
    [ -1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3],  # P
    [  1,-1, 1, 0,-1, 0,-1, 0,-1,-3,-3, 0,-2,-3,-1, 5, 2,-4,-2,-2],  # S
    [  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 2, 5,-3,-2, 0],  # T
    [ -3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1, 1,-4,-4,-3,15, 2,-3],  # W
    [ -2,-1,-2,-3,-3,-1,-2,-3, 2,-1,-1,-2, 0, 4,-3,-2,-2, 2, 8,-1],  # Y
    [  0,-3,-3,-4,-1,-3,-3,-4,-4, 4, 1,-3, 1,-1,-3,-2, 0,-3,-1, 5]   # V
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
  );

# The BLOSUM62 amino acid substitution matrix

@Constants::blosum62 = 
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
  ( [  4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],  # A
    [ -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],  # R
    [ -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],  # N
    [ -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],  # D
    [  0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],  # C
    [ -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],  # Q
    [ -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],  # E
    [  0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],  # G
    [ -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],  # H
    [ -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],  # I
    [ -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],  # L
    [ -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],  # K
    [ -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],  # M
    [ -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],  # F
    [ -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],  # P
    [  1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],  # S
    [  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],  # T
    [ -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],  # W
    [ -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],  # Y
    [  0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4]   # V
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
    );

# The amino acids, and a hashmap from amino acid to its index 0-19

my %aaIndex;

for (my $i=0; $i<20; $i++) {
  $aaIndex{$Constants::aas[$i]} = $aaIndex{lc($Constants::aas[$i])} = $i;
}

# The score of a pair of amino acids in a given matrix

sub score {
  my ($matrix, $aa1, $aa2) = @_;	# By ref, val, val
  return $$matrix[$aaIndex{$aa1}][$aaIndex{$aa2}];

}

# The maximum of a list of numbers

sub max {
  my $res = -1E300;
  foreach (@_) {
    if ($_ > $res) {
      $res = $_;
    }
  }
  return $res;
}

# Make random DNA sequence of specified length

sub randomdna {
  my ($len) = @_;	  # length of sequence to generate
  my $dna = "";
  foreach (1..$len) {
    my $p = rand(4);
    $dna = $dna . substr("ACGT", $p, 1);
  }
  return $dna;
}

1;   # Necessary for this module to load into scripts that use it.
