#! /usr/bin/perl 
#
# Description: Will loop through a given group of numbers and open
# their corresponding file in xcrysden using a save template file
#
# Input: 1) xcysden template file (save from a session)
#        2) File Root name. Must be of the form <Name><num>.xsf
#        2) a list of numbers corresponding to files
#
#------------------------------------------------------------------

use strict;
use warnings;


#Pull off the template file
my $file = shift @ARGV;

#file root
my $root = shift @ARGV;

#template file
open my $fh, '<', $file or die "ERROR: Cannot open $file $!";
my @file_lines = <$fh>;
close $fh;

#loop over the numbers
for my $num (@ARGV){

   #temp output file
   my $output = "output";
   open my $fout, '>', $output or die "ERROR: Cannot open $output $!";
   select $fout;

   for my $line (@file_lines){

      $line =~ s/(::scripting::open --xsf).*/$1 ${root}${num}.xsf/;

      print $line;

   }
   select STDOUT;
   close ($fout);

   system("xcrysden --script $output");
}
