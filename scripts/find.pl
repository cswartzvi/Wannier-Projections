#! /usr/bin/perl

use warnings;
use strict;

while ( my $line = <STDIN> ){
   my @temp = split /\s+/, $line;
   if ($temp[2] >= 0.10000) {
      print "$temp[1] $temp[2] \n";
   }
}
