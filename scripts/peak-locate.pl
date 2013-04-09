#! /usr/bin/perl
####################################################3
#
# Using the intKS RAW file to find the peaks for
# use in the 
#
#
####################################################3

use strict;
use warnings;
use 5.0120;

#intkS RAW File (This is a new file that has the eigenstate numbers of the x-axis);
my $file = shift @ARGV;

#--------------------------------------------------------
#Some major variables
my $sum  = 0.0;               #Sum of all Values
my $total_states = 0.0;       #Total Number of states in the intKS.dat RAW file
my $high_val = 0.0;           #Highest projection Peak
my $high_val_num = 0.0;       #State Number of Highest Projection Peak
my $peak_div;                 #Peak division threshold (What precent of max peak?) 

my @peaks;                    #Peak numbers
my @intKS_file;               #intKS RAW File (x - states, y - Raw Projections);

my $output = 'peaks.dat';     #main peak output
my $gnu = 'temp.gnu';         #gnuplot template file
my $temp = 'output.temp';     #output template file
#--------------------------------------------------------

open my $fh, '<', $file or die "ERROR: Cannot open $file $!";
@intKS_file = <$fh>;

#first loop through the file find the intKS_file find highest peak
for (@intKS_file){

   my ($num, $val) = split ' ', $_;

   if ($val > $high_val){
      $high_val = $val;
      $high_val_num = $num;
   }

   $total_states = $num if ($num > $total_states);

   $sum += $val;
   

}

print "Highest Value: $high_val at num: $high_val_num\n";
print "Highest State: $total_states\n";
print "Average Value: ", $high_val/$total_states, "\n";


#second loop through the file: Find all peaks
FIND_PEAKS: {
   print "\nFinding Peaks...\n";
   print "\nPlease enter the peak Threshold:";
   $peak_div = <STDIN>;
   print "\nPeak Threshold divsor : $peak_div\n";


   @peaks =();
   for (@intKS_file){

      
      my ($num, $val) = split ' ', $_;

      push @peaks, $num if ($val > (1/$peak_div)*($high_val));
   }

   #$Found Peaks
   print "Peaks: @peaks \n";

   open my $fout, '>', $temp;
   select $fout;
   for (@peaks){

      print " $_  ", 5*($high_val/$total_states), "\n";
   }
   select STDOUT;
   close $fout;

   #create gnuplot template
   open my $gnu_out, '>', $gnu or die "ERROR: Cannot Open $gnu $!";
   print $gnu_out "set term x11 title 'Peak division Threshold: $peak_div\n plot \'$file\' u 1:2 w lp, \'$temp\' u 1:2 w p lw 3\n pause -1";
   close $gnu_out;

   #create fork, open in gnuplot in a seperate terminal to check
   #if the peaks are correct 
   defined(my $pid = fork) or die "Cannot fork $!";
   unless ($pid){

      exec "gnome-terminal -x gnuplot $gnu";
      die "cannot exec date: $!";
   }
   waitpid($pid, 0);


   #prompt to make sure that the threshold is correct
   PROMPT: while (1){

      print "Are these peaks correct? [y/n]:";
      my $ans = <STDIN>;

      given ($ans) {
         when ( /y/i) { last PROMPT}
         when ( /n/i) { redo FIND_PEAKS}
         default      { next PROMPT}
      }
   }
}

#delete templates
unlink $temp, $gnu;

#Print the main peaks file
open my $main_out, '>', $output or die "ERROR: Cannot Open $output $!";
print $main_out "$_ \n" for (@peaks);
close $main_out;
