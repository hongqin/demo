#!/usr/bin/perl

# Needleman-Wunsch and Smith-Waterman alignment algorithms with
# graphic display of the dynamic programming matrix and traceback.
# Linear and affine gap costs.

# Peter Sestoft, Royal Veterinary and Agricultural University, Denmark
# Reference: http://www.itu.dk/people/sestoft/bsa.html
# sestoft@itu.dk * 2003-04-19, 2003-05-04, 2003-08-25, 2003-10-16

use strict;
use warnings;
use Constants;
use GD;
use CGI qw/:standard/;

my $author = "dina.kvl.dk";

# The sequences to align, the alignment type, the substitution matrix,
# indication whether the traceback should be shown, the gap opening
# cost and gap extension cost (the defaults are for testing)

my ($seq1, $seq2, $kind, $matrix, $showtraceback, $gapmodel, $d, $e) = 
    ("HEAGAWGHEE", "PAWHEAE", "global", 
     \@Constants::blosum50, 1, "linear", 8, 8);

if (param()) {
  $seq1 = param('seq1');
  $seq2 = param('seq2');
  $kind = param('kind');
  $matrix = 
      param('matrix') eq 'blosum45' ? \@Constants::blosum45 
    : param('matrix') eq 'blosum50' ? \@Constants::blosum50
    : param('matrix') eq 'blosum62' ? \@Constants::blosum62
    : \@Constants::blosum50;
  $showtraceback = param('showtraceback');
  ($gapmodel, $d, $e) = param('gapcost') =~ /([a-z]+)-(\d+)-(\d+)/;
}

# ----------------------------------------------------------------------
# Global constants

# The traceback matrices are indexed by (direction, row, column).

my @DIR = (1, 2, 3);
my $STOP = 0;

# Directions in the linear (2D) traceback: 
# 0=stop; 1=from North (above); 2=from Northwest; 3=from West (left)
my ($FROMN, $FROMNW, $FROMW) = @DIR;

# Directions in the affine (3D) traceback: 
my ($FROMM, $FROMIX, $FROMIY) = @DIR;

# Minus infinity

my $minusInf = -2111111111;     # Representable in 32 bits 2's compl.

# Color codes for the traceback
my ($RED, $BLUE, $GREEN) = (1, 2, 3);

# ----------------------------------------------------------------------
# Draw an image of the dynamic programming matrix and the traceback
# Linear gap costs
  
sub drawLinear {
  my ($seq1, $seq2, $Fref, $Bref, $seq1a, $seq2a) = @_;  
  my $font = GD::gdSmallFont;
  my ($fontw, $fonth) = ($font->width, $font->height);
  # Find widest number
  my $maxlen = 1;
  foreach my $rowref (@$Fref) {
    foreach my $x (@$rowref) {
      $maxlen = &max($maxlen, length($x));
    }
  }
  # Set up the dimensions and sizes of the image
  my $sep  = 10;
  my $marw = $fontw + $sep;
  my $marh = $fonth + $sep;
  my $colw = $maxlen * $fontw + $sep;
  my $rowh = $fonth + $sep;
  my $cols = length($seq1) + 1;
  my $rows = length($seq2) + 1;
  # Compute image dimensions, including matrix and alignment
  my $matrixw = $cols * $colw;
  my $matrixh = $rows * $rowh;
  my $alignw = $fontw * &max(length($seq1a), length($seq2a));
  my $imagew = &max($alignw, $marw + $matrixw + $sep);
  my $imageh = $sep + 2 * $fonth + $marh + $matrixh;
  # Allocate the image and the colors we need
  my $im = new GD::Image($imagew, $imageh);
  my $grey  = $im->colorAllocate(230,230,230);
  my $white = $im->colorAllocate(255,255,255);
  my $black = $im->colorAllocate(  0,  0,  0);
  my $red   = $im->colorAllocate(255,  0,  0);
  my $blue  = $im->colorAllocate(  0,  0,255);
  my $green = $im->colorAllocate(  0,255,  0);
  my @tracebackcolor = ($red, $blue, $green);  # Numbered 1, 2, 3
  # Write the sequences along the axes
  for (my $c=0; $c<$cols-1; $c++) {
    my $x = $marw + ($c+1) * $colw - $fontw/2 - 3;
    $im->string($font, $x, 2, substr($seq1, $c, 1), $black); 
  }
  for (my $r=0; $r<$rows-1; $r++) {
    my $y = $marh + $r * $rowh + $rowh/2;
    $im->string($font, 2, $y, substr($seq2, $r, 1), $black); 
  }
  { 
    # Draw borders around the matrix, and the internal dividing lines
    my $leftx  = $marw-$sep/2;
    my $rightx = $marw+$matrixw-5;
    my $topy   = $marh-$sep/2;
    my $boty   = $marh+$matrixh-5;
    $im->line(2,       $topy, $rightx, $topy, $black);
    $im->line(2,       $boty, $rightx, $boty, $black);
    $im->line($leftx,  2,     $leftx,  $boty, $black);
    $im->line($rightx, 2,     $rightx, $boty, $black);
    for (my $c=1; $c<$cols; $c++) {
      my $x = $marw + $c * $colw - 4;
      $im->filledRectangle($x-1, $topy+1, $x+1, $boty-1, $white); 
    }
    for (my $r=1; $r<$rows; $r++) {
      my $y = $marh + $r * $rowh - $sep/2;
      $im->filledRectangle($leftx+1, $y, $rightx-1, $y+2, $white); 
    }
  }
  # Write the elements of F and the arrows of B
  for (my $r=0; $r<$rows; $r++) {
    my $y0 = $marh + $r * $rowh;
    for (my $c=0; $c<$cols; $c++) {
      my $x0 = $marw + ($c + 1) * $colw - $sep;
      my $n = $$Fref[$r][$c];
      my $x = $x0 - length($n) * $fontw;
      $im->string($font, $x, $y0, $n, $black);
      if ($showtraceback == 1) {
	# Arrow pointing W (left)
	my $xa = $x0 - $colw;
	if ($$Bref[$FROMW][$r][$c]) {
	  my $color = $tracebackcolor[$$Bref[$FROMW][$r][$c] - 1];
	  $im->line($xa+1, $y0+$fonth/2+1, $xa+8, $y0+$fonth/2+1,  $color);
	  $im->line($xa+1, $y0+$fonth/2+1, $xa+3, $y0+$fonth/2-1, $color);
	  $im->line($xa+1, $y0+$fonth/2+1, $xa+3, $y0+$fonth/2+3, $color);
	}
	# Arrow pointing N (up)
	my $xb = $x0 - $sep;
	if ($$Bref[$FROMN][$r][$c]) {
	  my $color = $tracebackcolor[$$Bref[$FROMN][$r][$c] - 1];
	  $im->line($xb, $y0-1, $xb,   $y0-9, $color);
	  $im->line($xb, $y0-9, $xb-2, $y0-7, $color);
	  $im->line($xb, $y0-9, $xb+2, $y0-7, $color);
	}
	# Arrow pointing NW (up left)
	if ($$Bref[$FROMNW][$r][$c]) {
	  my $color = $tracebackcolor[$$Bref[$FROMNW][$r][$c] - 1];
	  $im->line($xa, $y0-10, $xa+10, $y0,    $color);
	  $im->line($xa, $y0-10, $xa,    $y0-7,  $color);
	  $im->line($xa, $y0-10, $xa+3,  $y0-10, $color);
	}
      }
    }
  }
  # Show the alignment
  my $alignx = ($imagew-$alignw)/2;
  $im->string($font, $alignx, $marh + $rows * $rowh + $sep/2,          
              $seq1a, $black);
  $im->string($font, $alignx, $marh + $rows * $rowh + $sep/2 + $fonth, 
              $seq2a, $black);
  # Put the author in there
  my $tinyfont = GD::gdTinyFont;
  my ($tinyfontw, $tinyfonth) = ($tinyfont->width, $tinyfont->height);
  my $authory = ($matrixh+$tinyfontw*length($author))/2 + $marh;
  $im->stringUp($tinyfont, $imagew-$tinyfonth, $authory, $author, $black);
  print $im->png;
}

# ----------------------------------------------------------------------
# The Needleman-Wunsch global alignment algorithm, linear gap costs
# Input: The sequences $x and $y to align
# Output: references to the F and B matrices, and the aligned sequences

sub globalAlignLinear {
  my ($matrix, $x, $y) = @_;            # By ref, val, val
  my ($n, $m) = (length($x), length($y));
  # The dynamic programming matrix
  my @F; 
  for (my $j=0; $j<=$m; $j++) {
    $F[$j] = [(0) x ($n+1)];
  }
  # The traceback matrix
  my @B; 
  foreach my $dir (@DIR) {
    for (my $j=0; $j<=$m; $j++) {
      $B[$dir][$j] = [(0) x ($n+1)];
    }
  }
  # Initialize upper and left-hand borders of F and B matrices
  for (my $i=1; $i<=$n; $i++) {
    $F[0][$i] = -$e * $i;
    $B[$FROMW][0][$i] = $RED;
  }
  for (my $j=1; $j<=$m; $j++) {
    $F[$j][0] = -$e * $j;
    $B[$FROMN][$j][0] = $RED;
  }
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      my $s = &score($matrix, substr($x, $i-1, 1), substr($y, $j-1, 1));
      my $val = &max($F[$j-1][$i-1]+$s, 
                     $F[$j][$i-1]-$e, 
                     $F[$j-1][$i]-$e);
      $F[$j][$i] = $val;
      # Record all traceback directions
      if ($val == $F[$j-1][$i-1]+$s) {
        $B[$FROMNW][$j][$i] = $RED;
      } 
      if ($val == $F[$j][$i-1]-$e) {
        $B[$FROMW][$j][$i] = $RED;
      } 
      if ($val == $F[$j-1][$i]-$e) {
        $B[$FROMN][$j][$i] = $RED;
      } 
    }
  }
  &markReachable2(\@B, $n, $m);
  return (\@F, \@B, &traceback2($x, $y, \@B, $n, $m));
}

# ----------------------------------------------------------------------
# The Smith-Waterman local alignment algorithm, linear gap costs
# Input: The sequences $x and $y to align
# Output: references to the F and B matrices, and the aligned sequences

sub localAlignLinear {
  my ($matrix, $x, $y) = @_;            # By ref, val, val
  my ($n, $m) = (length($x), length($y));
  # The dynamic programming matrix; also correctly initializes borders
  my @F; 
  for (my $j=0; $j<=$m; $j++) {
    $F[$j] = [(0) x ($n+1)];
  }
  # The traceback matrix; also correctly initializes borders
  my @B; 
  foreach my $dir (@DIR) {
    for (my $j=0; $j<=$m; $j++) {
      $B[$dir][$j] = [($STOP) x ($n+1)];
    }
  }
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      my $s = &score($matrix, substr($x, $i-1, 1), substr($y, $j-1, 1));
      my $val = &max(0, 
                     $F[$j-1][$i-1]+$s, 
                     $F[$j][$i-1]-$e, 
                     $F[$j-1][$i]-$e);
      $F[$j][$i] = $val;
      # Record all traceback directions (none if we restart at score 0):
      if ($val == $F[$j-1][$i-1]+$s) {
        $B[$FROMNW][$j][$i] = $RED;
      } 
      if ($val == $F[$j][$i-1]-$e) {
        $B[$FROMW][$j][$i] = $RED;
      } 
      if ($val == $F[$j-1][$i]-$e) {
        $B[$FROMN][$j][$i] = $RED;
      } 
    }
  }
  # Find maximal score in matrix
  my $vmax = 0;
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      $vmax = &max($vmax, $F[$j][$i]);
    }
  }  
  my ($jmax, $imax) = (0, 0);
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      if ($vmax == $F[$j][$i]) {
        &markReachable2(\@B, $i, $j);
        $jmax = $j;
        $imax = $i;
      }
    }
  }  
  return (\@F, \@B, &traceback2($x, $y, \@B, $imax, $jmax));
}

# ----------------------------------------------------------------------
# Common subroutines for linear gap cost routines

# Reconstruct the alignment from the traceback, backwards, from ($i, $j)

sub traceback2 {
  my ($x, $y, $B, $i, $j) = @_;         # B by reference
  my ($xAlign, $yAlign) = ("", "");
  while ($$B[$FROMN][$j][$i] || $$B[$FROMW][$j][$i] || $$B[$FROMNW][$j][$i]) {
    if ($$B[$FROMN][$j][$i]) {
      $$B[$FROMN][$j][$i] = $GREEN;
      $xAlign .= "-"; 
      $yAlign .= substr($y, $j-1, 1);
      $j--;
    } elsif ($$B[$FROMW][$j][$i]) {
      $$B[$FROMW][$j][$i] = $GREEN;
      $xAlign .= substr($x, $i-1, 1);
      $yAlign .= "-"; 
      $i--;
    } elsif ($$B[$FROMNW][$j][$i]) {
      $$B[$FROMNW][$j][$i] = $GREEN;
      $xAlign .= substr($x, $i-1, 1);
      $yAlign .= substr($y, $j-1, 1);
      $i--; $j--;
    }
  }
  # Warning: these expressions cannot be inlined in the list
  $xAlign = reverse $xAlign;
  $yAlign = reverse $yAlign;
  return ($xAlign, $yAlign);
}

# Mark all traceback arrows reachable from a ($i, $j)

sub markReachable2 {
  my ($B, $i, $j) = @_;         # B by reference
  if ($$B[$FROMN][$j][$i] == $RED) {
    $$B[$FROMN][$j][$i] = $BLUE;
    &markReachable2($B, $i, $j-1);
  } 
  if ($$B[$FROMW][$j][$i] == $RED) {
    $$B[$FROMW][$j][$i] = $BLUE;
    &markReachable2($B, $i-1, $j);
  } 
  if ($$B[$FROMNW][$j][$i] == $RED) {
    $$B[$FROMNW][$j][$i] = $BLUE;
    &markReachable2($B, $i-1, $j-1);
  }
}

# ----------------------------------------------------------------------
# Replace the number $minusInf by the string "-Inf" in a matrix

sub replaceMinusInf {
  my ($matref) = @_;            # By reference
  foreach my $rowref (@$matref) {
    foreach my $x (@$rowref) {
      if ($x == $minusInf) {
        $x = "-Inf";
      }
    }
  }
}

# ----------------------------------------------------------------------
# Draw an image of the dynamic programming matrix and the traceback
# Affine gap costs
  
sub drawAffine {
  my ($seq1, $seq2, $Mref, $Ixref, $Iyref, $Bref, $seq1a, $seq2a) = @_;  
  my $font = GD::gdSmallFont;
  my ($fontw, $fonth) = ($font->width, $font->height);
  my $cols = length($seq1) + 1;
  my $rows = length($seq2) + 1;
  # Replace $minusInf by "-Inf"
  &replaceMinusInf($Mref);
  &replaceMinusInf($Ixref);
  &replaceMinusInf($Iyref);
  # Find widest number, taking into account the displacement of M and Ix
  my $maxlen = 1;
  for (my $r=0; $r<$rows; $r++) {
    for (my $c=0; $c<$cols; $c++) {
      $maxlen = &max($maxlen, 
                     length($$Iyref[$r][$c]),
                     1+length($$Mref[$r][$c]),
                     2+length($$Ixref[$r][$c]));
    }
  }
  # Set up the dimensions and sizes of the image
  my $sep  = 10;
  my $marw = $fontw + $sep;
  my $marh = $fonth + $sep;
  my $colw = $maxlen * $fontw + $sep;
  my $rowh = 3 * $fonth + $sep;
  # Compute image dimensions, including matrix and alignment
  my $matrixw = $cols * $colw;
  my $matrixh = $rows * $rowh;
  my $alignw = $fontw * &max(length($seq1a), length($seq2a));
  my $imagew = &max($alignw, $marw + $matrixw + $sep);
  my $imageh = $sep + 2 * $fonth + $marh + $matrixh;
  # Allocate the image and the colors we need
  my $im = new GD::Image($imagew, $imageh);
  my $grey  = $im->colorAllocate(230,230,230);
  my $white = $im->colorAllocate(255,255,255);
  my $black = $im->colorAllocate(  0,  0,  0);
  my $red   = $im->colorAllocate(255,  0,  0);
  my $blue  = $im->colorAllocate(  0,  0,255);
  my $green = $im->colorAllocate(  0,255,  0);
  my @tracebackcolor = ($red, $blue, $green);  # Numbered 1, 2, 3
  # Write the sequences along the axes
  for (my $c=1; $c<$cols; $c++) {
    my $x = $marw + $c * $colw - $fontw/2 - 3;
    $im->string($font, $x, 2, substr($seq1, $c-1, 1), $black); 
  }
  for (my $r=1; $r<$rows; $r++) {
    my $y = $marh + $r * $rowh - $fonth + 2;
    $im->string($font, 2, $y, substr($seq2, $r-1, 1), $black); 
  }
  { 
    # Draw borders around the matrix, and the internal dividing lines
    my $leftx  = $marw-$sep/2;
    my $rightx = $marw+$matrixw-3;
    my $topy   = $marh-$sep/2;
    my $boty   = $marh+$matrixh;
    $im->line(2,       $topy, $rightx, $topy, $black);
    $im->line(2,       $boty, $rightx, $boty,  $black);
    $im->line($leftx,  2,     $leftx,  $boty, $black);
    $im->line($rightx, 2,     $rightx, $boty,  $black);
    for (my $c=1; $c<$cols; $c++) {
      my $x = $marw + $c * $colw - 4;
      $im->filledRectangle($x-1, $topy+1, $x+1, $boty-1, $white); 
    }
    for (my $r=1; $r<$rows; $r++) {
      my $y = $marh + $r * $rowh - $sep/2;
      $im->filledRectangle($leftx+1, $y, $rightx-1, $y+2, $white); 
    }
  }
  # Write the elements of M, Ix, Iy and the arrows of B
  for (my $r=0; $r<$rows; $r++) {
    # Top coordinate of this row
    my $y0 = $marh + $r * $rowh;
    for (my $c=0; $c<$cols; $c++) {
      # Right-hand coordinate of this column
      my $x0 = $marw + ($c + 1) * $colw - $sep;
      # Show the Iy, M and Ix matrix entries.  Cell layout:
      #     Iy
      #     M
      #   Ix
      my $n = $$Iyref[$r][$c];
      my $x = $x0 - length($n) * $fontw;
      $im->string($font, $x, $y0, $n, $black);
      $n = $$Mref[$r][$c];
      $x = $x0 - (1 + length($n)) * $fontw;
      $im->string($font, $x, $y0+$fonth, $n, $black);
      $n = $$Ixref[$r][$c];
      $x = $x0 - (2 + length($n)) * $fontw;
      $im->string($font, $x, $y0+2*$fonth, $n, $black);
      if ($showtraceback == 1) {
	{ # Arrows pointing up from Iy 
	  my $taily = $y0;
	  # Arrow pointing to Iy
	  if ($$Bref[$FROMIY][3][$r][$c]) {
	    my $color = $tracebackcolor[$$Bref[$FROMIY][3][$r][$c] - 1];
	    my $headx = $x0-0.5*$fontw-1;
	    my $heady = $y0-$sep-2*$fonth+4;
	    $im->line($headx, $heady, $headx,   $taily,   $color);
	    $im->line($headx, $heady, $headx-3, $heady+3, $color);
	    $im->line($headx, $heady, $headx+3, $heady+3, $color);
	  }
	  # Arrow pointing to M
	  if ($$Bref[$FROMM][3][$r][$c]) {
	    my $color = $tracebackcolor[$$Bref[$FROMM][3][$r][$c] - 1];
	    my $headx = $x0-1.5*$fontw-1;
	    my $heady = $y0-$sep-$fonth-1;
	    $im->line($headx, $heady, $headx,   $taily,   $color);
	    $im->line($headx, $heady, $headx-3, $heady+3, $color);
	    $im->line($headx, $heady, $headx+3, $heady+3, $color);
	  }
	  # Arrow pointing to Ix
	  if ($$Bref[$FROMIX][3][$r][$c]) {
	    my $color = $tracebackcolor[$$Bref[$FROMIX][3][$r][$c] - 1];
	    my $headx = $x0-2.5*$fontw-1;
	    my $heady = $y0-$sep-1;
	    $im->line($headx, $heady, $headx,   $taily,   $color);
	    $im->line($headx, $heady, $headx-3, $heady+3, $color);
	    $im->line($headx, $heady, $headx+3, $heady+3, $color);
	  }
	}
	{ # Arrows pointing up left from M 
	  my $headx = $x0 - $colw;
	  my $tailx = $headx+$fontw+13;
	  my $taily = $y0+1.5*$fonth-5;
	  # Arrow pointing to Iy
	  if ($$Bref[$FROMIY][1][$r][$c]) {
	    my $heady = $y0-$rowh+$fonth/2+6;
	    my $color = $tracebackcolor[$$Bref[$FROMIY][1][$r][$c] - 1];
	    $im->line($headx, $heady, $tailx,   $taily,   $color);
	    $im->line($headx, $heady, $headx-2, $heady+4, $color);
	    $im->line($headx, $heady, $headx+4, $heady,   $color);
	  }
	  # Arrow pointing to M
	  if ($$Bref[$FROMM][1][$r][$c]) {
	    my $color = $tracebackcolor[$$Bref[$FROMM][1][$r][$c] - 1];
	    my $heady = $y0-$rowh+1.5*$fonth+7;
	    my $hx = $headx-1;
	    $im->line($hx, $heady, $tailx, $taily,   $color);
	    $im->line($hx, $heady, $hx,    $heady+4, $color);
	    $im->line($hx, $heady, $hx+4,  $heady,   $color);
	  }
	  # Arrow pointing to Ix
	  if ($$Bref[$FROMIX][1][$r][$c]) {
	    my $color = $tracebackcolor[$$Bref[$FROMIX][1][$r][$c] - 1];
	    my $hx = $headx;
	    my $heady = $y0-$rowh+2.5*$fonth+5;
	    $im->line($hx, $heady, $tailx, $taily,   $color);
	    $im->line($hx, $heady, $hx,    $heady+4, $color);
	    $im->line($hx, $heady, $hx+4,  $heady,   $color);
	  }
	}
	{ # Arrows pointing left from Ix 
	  my $headx = $x0 - $colw;
	  my $tailx = $headx+12;
	  my $taily = $y0+2.5*$fonth+1;
	  # Arrow pointing to Iy
	  if ($$Bref[$FROMIY][2][$r][$c]) {
	    my $heady = $y0+$fonth/2+6;
	    my $color = $tracebackcolor[$$Bref[$FROMIY][2][$r][$c] - 1];
	    $im->line($headx, $heady, $tailx,   $taily,   $color);
	    $im->line($headx, $heady, $headx-2, $heady+4, $color);
	    $im->line($headx, $heady, $headx+4, $heady,   $color);
	  }
	  # Arrow pointing to M
	  if ($$Bref[$FROMM][2][$r][$c]) {
	    my $color = $tracebackcolor[$$Bref[$FROMM][2][$r][$c] - 1];
	    my $heady = $y0+1.5*$fonth+7;
	    my $hx = $headx-1;
	    $im->line($hx, $heady, $tailx, $taily,   $color);
	    $im->line($hx, $heady, $hx,    $heady+4, $color);
	    $im->line($hx, $heady, $hx+4,  $heady,   $color);
	  }
	  # Arrow pointing to Ix
	  if ($$Bref[$FROMIX][2][$r][$c]) {
	    my $color = $tracebackcolor[$$Bref[$FROMIX][2][$r][$c] - 1];
	    my $hx = $headx-1.5*$fontw;
	    my $heady = $y0+2.5*$fonth+1;
	    $im->line($hx, $heady, $tailx, $heady,   $color);
	    $im->line($hx, $heady, $hx+3,  $heady-3, $color);
	    $im->line($hx, $heady, $hx+3,  $heady+3, $color);
	  }
	}
      }
    }
  }
  # Show the alignment
  my $alignx = ($imagew-$alignw)/2;
  $im->string($font, $alignx, $marh + $rows * $rowh + $sep/2,          
              $seq1a, $black);
  $im->string($font, $alignx, $marh + $rows * $rowh + $sep/2 + $fonth, 
              $seq2a, $black);
  # Put the author in there
  my $tinyfont = GD::gdTinyFont;
  my ($tinyfontw, $tinyfonth) = ($tinyfont->width, $tinyfont->height);
  my $authory = ($matrixh+$tinyfontw*length($author))/2 + $marh;
  $im->stringUp($tinyfont, $imagew-$tinyfonth, $authory, $author, $black);
  print $im->png;
}


# ----------------------------------------------------------------------
# The Needleman-Wunsch global alignment algorithm, affine gap costs
# Input: The sequences $x and $y to align
# Output: references to the matrices M, Ix, Iy, B, and the aligned sequences

sub globalAlignAffine {
  my ($matrix, $x, $y) = @_;
  my ($n, $m) = (length($x), length($y));
  # Initialize upper and left-hand borders
  # M represent an aa/aa match; 
  # Ix represents insertions in x (gaps in y); 
  # Iy represents insertions in y (gaps in x); 
  # The traceback now points to the matrix (M, Ix, Iy) from which the
  # maximum was obtained: $FROMM=1, $FROMIX=2, $FROMIY=3
  # B[$dir][1] is the traceback for M; 
  # B[$dir][2] is the traceback for Ix; 
  # B[$dir][3] is the traceback for Iy
  my (@M, @Ix, @Iy, @B);
  $M[0][0] = 0;
  $Ix[0][0] = $Iy[0][0] = $minusInf;
  foreach my $dir (@DIR) {
    for (my $j=0; $j<=$m; $j++) {
      for (my $k=1; $k<=3; $k++) {
        $B[$dir][$k][$j] = [($STOP) x ($n+1)];
      }
    }
  }
  for (my $i=1; $i<=$n; $i++) {
    $Ix[0][$i] = - $d - $e * ($i-1);
    $B[$FROMIX][2][0][$i] = $RED;
    $Iy[0][$i] = $minusInf;
    $M[0][$i] = $minusInf;
  }
  for (my $j=1; $j<=$m; $j++) {
    $Iy[$j][0] = - $d - $e * ($j-1);
    $B[$FROMIY][3][$j][0] = $RED;
    $Ix[$j][0] = $minusInf;
    $M[$j][0] = $minusInf;
  }
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      my $s = &score($matrix, substr($x, $i-1, 1), substr($y, $j-1, 1));
      my $val = &max($M[$j-1][$i-1]+$s, 
                     $Ix[$j-1][$i-1]+$s, 
                     $Iy[$j-1][$i-1]+$s);
      $M[$j][$i] = $val;
      if ($val == $M[$j-1][$i-1]+$s) {
        $B[$FROMM][1][$j][$i] = $RED; 
      } 
      if ($val == $Ix[$j-1][$i-1]+$s) {
        $B[$FROMIX][1][$j][$i] = $RED; 
      } 
      if ($val == $Iy[$j-1][$i-1]+$s) {
        $B[$FROMIY][1][$j][$i] = $RED; 
      } 
      $val = &max($M[$j][$i-1]-$d, $Ix[$j][$i-1]-$e, $Iy[$j][$i-1]-$d);  
      $Ix[$j][$i] = $val;
      if ($val == $M[$j][$i-1]-$d) {
        $B[$FROMM][2][$j][$i] = $RED;
      } 
      if ($val == $Ix[$j][$i-1]-$e) {
        $B[$FROMIX][2][$j][$i] = $RED;
      } 
      if ($val == $Iy[$j][$i-1]-$d) {
        $B[$FROMIY][2][$j][$i] = $RED;
      }      
      $val = &max($M[$j-1][$i]-$d, $Iy[$j-1][$i]-$e, $Ix[$j-1][$i]-$d);  
      $Iy[$j][$i] = $val;
      if ($val == $M[$j-1][$i]-$d) {
        $B[$FROMM][3][$j][$i] = $RED;
      } 
      if ($val == $Iy[$j-1][$i]-$e) {
        $B[$FROMIY][3][$j][$i] = $RED;
      } 
      if ($val == $Ix[$j-1][$i]-$d) {
        $B[$FROMIX][3][$j][$i] = $RED;
      }      
    }
  }
  # Find the matrix (@M, @Ix or @Iy) whose ($m,$n) cell has highest score:
  my ($kmax, $vmax) = (1, $M[$m][$n]);
  if ($vmax < $Ix[$m][$n]) {
    $vmax = $Ix[$m][$n];
    $kmax = 2;
  }
  if ($vmax < $Iy[$m][$n]) {
    $vmax = $Iy[$m][$n];
    $kmax = 3;
  }
  if ($M[$m][$n] == $vmax) {
    &markReachable3(\@B, 1, $n, $m);
  }
  if ($Ix[$m][$n] == $vmax) {
    &markReachable3(\@B, 2, $n, $m);
  }
  if ($Iy[$m][$n] == $vmax) {
    &markReachable3(\@B, 3, $n, $m);
  }
  return (\@M, \@Ix, \@Iy, \@B, &traceback3($x, $y, \@B, $kmax, $n, $m));
}


# ----------------------------------------------------------------------
# The Smith-Waterman local alignment algorithm, affine gap costs
# Input: The sequences $x and $y to align
# Output: references to the matrices M, Ix, Iy, B, and the aligned sequences

sub localAlignAffine {
  my ($matrix, $x, $y) = @_;
  my ($n, $m) = (length($x), length($y));
  # Initialize upper and left-hand borders
  # M represent an aa/aa match; 
  # Ix represents insertions in x (gaps in y); 
  # Iy represents insertions in y (gaps in x); 
  # The traceback now points to the matrix (M, Ix, Iy) from which the
  # maximum was obtained: $FROMM=1, $FROMIX=2, $FROMIY=3
  # B[$dir][1] is the traceback for M; 
  # B[$dir][2] is the traceback for Ix; 
  # B[$dir][3] is the traceback for Iy
  my (@M, @Ix, @Iy, @B);
  for (my $j=0; $j<=$m; $j++) {
    $M[$j] = [(0) x ($n+1)];
    $Ix[$j] = [($minusInf) x ($n+1)];
    $Iy[$j] = [($minusInf) x ($n+1)];
  }
  # The traceback matrix; also correctly initializes borders
  foreach my $dir (@DIR) {
    for (my $j=0; $j<=$m; $j++) {
      for (my $k=1; $k<=3; $k++) {
        $B[$dir][$k][$j] = [($STOP) x ($n+1)];
      }
    }
  }
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      my $s = &score($matrix, substr($x, $i-1, 1), substr($y, $j-1, 1));
      my $val = &max(0, 
                     $M[$j-1][$i-1]+$s, 
                     $Ix[$j-1][$i-1]+$s, 
                     $Iy[$j-1][$i-1]+$s);
      $M[$j][$i] = $val;
      if ($val == $M[$j-1][$i-1]+$s) {
        $B[$FROMM][1][$j][$i] = $RED; 
      } 
      if ($val == $Ix[$j-1][$i-1]+$s) {
        $B[$FROMIX][1][$j][$i] = $RED; 
      } 
      if ($val == $Iy[$j-1][$i-1]+$s) {
        $B[$FROMIY][1][$j][$i] = $RED; 
      } 
      $val = &max($M[$j][$i-1]-$d, $Ix[$j][$i-1]-$e, $Iy[$j][$i-1]-$d);  
      $Ix[$j][$i] = $val;
      if ($val == $M[$j][$i-1]-$d) {
        $B[$FROMM][2][$j][$i] = $RED;
      } 
      if ($val == $Ix[$j][$i-1]-$e) {
        $B[$FROMIX][2][$j][$i] = $RED;
      } 
      if ($val == $Iy[$j][$i-1]-$d) {
        $B[$FROMIY][2][$j][$i] = $RED;
      }      
      $val = &max($M[$j-1][$i]-$d, $Iy[$j-1][$i]-$e, $Ix[$j-1][$i]-$d);  
      $Iy[$j][$i] = $val;
      if ($val == $M[$j-1][$i]-$d) {
        $B[$FROMM][3][$j][$i] = $RED;
      } 
      if ($val == $Iy[$j-1][$i]-$e) {
        $B[$FROMIY][3][$j][$i] = $RED;
      } 
      if ($val == $Ix[$j-1][$i]-$d) {
        $B[$FROMIX][3][$j][$i] = $RED;
      }      
    }
  }
  # Find maximal score in matrices
  my $vmax = 0;
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      $vmax = &max($vmax, $M[$j][$i], $Ix[$j][$i], $Iy[$j][$i]);
    }
  }  
  my ($kmax, $jmax, $imax) = (0, 0);
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      if ($vmax == $M[$j][$i]) {
        &markReachable3(\@B, 1, $i, $j);
        $kmax = 1;
        $jmax = $j;
        $imax = $i;
      }
      if ($vmax == $Ix[$j][$i]) {
        &markReachable3(\@B, 2, $i, $j);
        $kmax = 2;
        $jmax = $j;
        $imax = $i;
      }
      if ($vmax == $Iy[$j][$i]) {
        &markReachable3(\@B, 3, $i, $j);
        $kmax = 3;
        $jmax = $j;
        $imax = $i;
      }
    }
  }  
  return (\@M, \@Ix, \@Iy, \@B, &traceback3($x, $y, \@B, $kmax, $imax, $jmax));
}


# ------------------------------------------------------------
# Common subroutines for affine gap cost alignment
# Reconstruct the alignment from the traceback, backwards, 
# and mark green the path actually taken

sub traceback3 {
  my ($x, $y, $B, $k, $i, $j) = @_;   # B by reference
  my ($xAlign, $yAlign) = ("", "");
  while ($$B[$FROMM][$k][$j][$i] != 0 
         || $$B[$FROMIX][$k][$j][$i] != 0 
         || $$B[$FROMIY][$k][$j][$i] != 0) {
    my $nextk;
    # Colour green the path that was actually taken
    if ($$B[$FROMIY][$k][$j][$i]) {
      $$B[$FROMIY][$k][$j][$i] = $GREEN;
      $nextk = 3;       # From Iy
    } elsif ($$B[$FROMIX][$k][$j][$i]) {
      $$B[$FROMIX][$k][$j][$i] = $GREEN;
      $nextk = 2;       # From Ix
    } elsif ($$B[$FROMM][$k][$j][$i]) {
      $$B[$FROMM][$k][$j][$i] = $GREEN;
      $nextk = 1;       # From M
    } 
    if ($k == 1) {              # We're in the M matrix
      $xAlign .= substr($x, $i-1, 1);
      $yAlign .= substr($y, $j-1, 1);
      $i--; $j--;
    } elsif ($k == 2) {         # We're in the Ix matrix
      $xAlign .= substr($x, $i-1, 1);
      $yAlign .= "-"; 
      $i--;
    } elsif ($k == 3) {         # We're in the Iy matrix
      $xAlign .= "-"; 
      $yAlign .= substr($y, $j-1, 1);
      $j--;
    }       
    $k = $nextk;
  }
  # Warning: these expressions cannot be inlined in the list
  $xAlign = reverse $xAlign;
  $yAlign = reverse $yAlign;
  return ($xAlign, $yAlign);
}


# Mark blue all (affine) traceback arrows reachable from ($k, $i, $j)

sub markReachable3 {
  my ($B, $k, $i, $j) = @_;             # B by reference
  foreach my $dir (@DIR) {
    if ($$B[$dir][$k][$j][$i] == $RED) {
      $$B[$dir][$k][$j][$i] = $BLUE;
      if ($k == 1) {                    # We're in the M matrix
        &markReachable3($B, $dir, $i-1, $j-1);
      } elsif ($k == 2) {               # We're in the Ix matrix
        &markReachable3($B, $dir, $i-1, $j);
      } elsif ($k == 3) {               # We're in the Iy matrix
        &markReachable3($B, $dir, $i,   $j-1);
      }
    }
  }
}

#----------------------------------------------------------------------
# For debugging only
# Print a given matrix (array of array of number)
# The matrix must be passed by reference as in printmatrix(\@blosum50)

sub printmatrix {
  my ($title, $matrix) = @_;
  print "-" x 70, "\n$title:\n";
  foreach my $row (@$matrix) {
    foreach my $score (@$row) {
      if ($score == $minusInf) {
        $score = "-Inf";
        printf("%4s", $score);
      } else {
        printf("%4d", $score);
      }
    }
    print "\n";
  }
  print "-" x 70, "\n";
}

# sub printmatrix { }

# ----------------------------------------------------------------------
# For input sanity checking

sub badaa {
  my ($seq) = @_;
  $seq =~ tr/ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv//d;
  return $seq;
}

sub onlynucl {
  my ($seq) = @_;
  my $nuc = $seq =~ tr/ACGTUacgtu//;
  return $nuc > 12 && $nuc == length($seq);
}

sub errorpage {
  my ($msg) = @_;
  print header(-type=>'text/html');
  print "<HTML><BODY>\n
         <H2>Input error:</H2>\n
         $msg\n
         <P>\n
         Go back and try again!\n
         </BODY></HTML>\n";
}

# ----------------------------------------------------------------------
# Main program: 
# Check arguments and call the alignment and graph drawing functions

my $bad1 = &badaa($seq1);
my $bad2 = &badaa($seq2);

my $cellmax = 10000;

if ($bad1 ne "" || $bad2 ne "") {
  &errorpage("Your input sequences contain illegal characters: '$bad1$bad2'.
              <BR>They should contain amino acids only: ARNDCQEGHILKMFPSTWYV");
} elsif (&onlynucl($seq1 . $seq2)) {
  &errorpage("It is not meaningful to align nucleotide sequences (DNA, RNA)
              using this program.  It is for amino acids: ARNDCQEGHILKMFPSTWYV");
} elsif (length($seq1) * length($seq2) > $cellmax) {
  &errorpage("Your sequences are too long for this server.  The product of 
              the sequence lengths must be at most $cellmax.");
} else {
  print header(-type=>'image/png');
  if ($gapmodel eq "linear") {
    if ($kind eq "global") { 
      my ($F, $B, $xa, $ya) = &globalAlignLinear($matrix, $seq1, $seq2);
      &drawLinear($seq1, $seq2, $F, $B, $xa, $ya);
    } elsif ($kind eq "local") { 
      my ($F, $B, $xa, $ya) = &localAlignLinear($matrix, $seq1, $seq2);
      &drawLinear($seq1, $seq2, $F, $B, $xa, $ya);
    }
  } elsif  ($gapmodel eq "affine") {
    if ($kind eq "global") { 
      my ($M, $Ix, $Iy, $B, $xa, $ya) 
	  = &globalAlignAffine($matrix, $seq1, $seq2);
      &drawAffine($seq1, $seq2, $M, $Ix, $Iy, $B, $xa, $ya);
    } elsif ($kind eq "local") { 
      my ($M, $Ix, $Iy, $B, $xa, $ya)  
	  = &localAlignAffine($matrix, $seq1, $seq2);
      &drawAffine($seq1, $seq2, $M, $Ix, $Iy, $B, $xa, $ya);
    }
  }
}
