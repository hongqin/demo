
my $infile1 = 'marker.txt';  my $outfile1 = 'marker2.txt'; 
my $infile2 = 'trait.txt';   my $outfile2 = 'trait2.txt'; 

open (IN, "<$infile1"); open(OUT, ">$outfile1");
while( my $line=<IN> ) {
 $line =~ s/NA/-/g;
 print OUT $line; 
}
close(IN); close (OUT); 

my $debug = 10; 

open (IN, "<$infile2"); open(OUT, ">$outfile2");
while( my $line=<IN> ) {
 #$line =~ s/NA/-/g;
 chomp $line; 
 my @els = split( /\t/, $line ); 
 for(my $i=0; $i<=$#els; $i=$i+1 ) {
   if (defined $els[$i] ) {
     if ($els[$i] =~ /^\s*$/) { $els[$i]='NA'; }
     if ($i==$#els) {print "$els[0]-->$els[$i]\n"}
   } else {
       $els[$i] = 'NA'; 
   }
   if ($i ==$#els ) {
       print OUT $els[$i]."\n"; 
   } else {
       print OUT $els[$i]."\t"; 
   }
 }

}
close(IN); close (OUT); 

