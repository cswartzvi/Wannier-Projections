#! /usr/bin/perl
###########################################################
#
# Input
# 1) KS Files root
# 2) KS file ext
# 3) peaks file
#
###########################################################

use strict;
use warnings;

#root and ext of the KS Files
my $root = shift @ARGV;
my $ext = shift @ARGV;

my $tar_line = "tar -cvzf projections.tgz eig.dat intKS.* peaks.dat ";
for (<>){

   my $num = (split)[1];

   $tar_line .= " ".$root.$num.$ext." ";

}

system ("$tar_line");
