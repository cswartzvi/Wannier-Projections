#! /usr/bin/perl
#

use strict;
use warnings;

#first input is the gnuplot template file
my $temp_file = shift @ARGV;
open my $fh, '<', $temp_file or die "ERROR: $temp_file does not exist!! $!";
my @temp = <$fh>;
close ($fh);


my $output = 'plot-temp';

#The rst of the arguments are for the number of files
for my $file (glob 'output.dat*' ){
   
   (my $num = $file) =~ s/output.dat//;

   print "File: $file, Number: $num\n";

   open my $fout, '>', $output or die "ERROR: $output does not exist!! $!";
   select $fout;

   for my $line (@temp){
      
      $line =~ s/.dat\d+/.dat$num/;
      $line =~ s/(plot-)\d+/$1$num/;

      print $line;
   }
   select STDOUT;
   close ($fout);

   system("gnuplot -e \" load \'$output\' \" ");
}

