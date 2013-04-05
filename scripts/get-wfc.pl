#! /usr/bin/perl

use strict;
use warnings;

#Get the invocation arugments
my ($num, $pos_file, $wan_file, $alat) = @ARGV;

open my $fh, '<', "$pos_file";
my @pos_file = <$fh>;
my @ion_pos = split ' ', $pos_file[$num-1];
close($fh);


my @wan_dist;
open $fh, '<', "$wan_file";

my $wan_num;
while (my $line = <$fh>){ 

   $wan_num++;

   my @wan_pos = split ' ', $line;

   my @dr;
   for my $i (0 .. 2){
      $dr[$i] = $wan_pos[$i] - $ion_pos[$i];
      $dr[$i] = $dr[$i] - &nint($dr[$i]/$alat)*$alat;
   }

   my $r = sqrt($dr[0]**2 + $dr[1]**2 + $dr[2]**2);
   
   #if ( $r < 1.5){
   push @wan_dist, sprintf( "%15.10f %d",$r,$wan_num);
   #}
      
} 
close($fh);

#print "$_ $wan_dist[$_]\n" for (0 .. $#wan_dist);

my @sorted_dist = sort {
   my $anum = (split ' ', $a)[0];
   my $bnum = (split ' ', $b)[0];

   $anum <=> $bnum;
} @wan_dist;

print "$_ \n" for (@sorted_dist);
      
sub nint {  my $x = $_[0]; 
   my $n = int($x);
   if ( $x > 0 ) {    
      if ( $x-$n > 0.5) {      
         return $n+1;
      }
      else {
         return $n;
      }
   }
   else {    
      if ( $n-$x > 0.5) {      
         return $n-1;
      }
      else {      
         return $n;
      }
   } 
}
   
