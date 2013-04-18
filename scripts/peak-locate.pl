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

#intkS RAW File (This is a new file that has the eigenstate numbers of the x-axis);
my $file_raw = shift @ARGV;
my $file = shift @ARGV;

#--------------------------------------------------------
#Some major variables
my $sum  = 0.0;               #Sum of all Values
my $total_states = 0.0;       #Total Number of states in the intKS.dat RAW file
my $high_val = 0.0;           #Highest projection Peak
my $high_val_num = 0.0;       #State Number of Highest Projection Peak
my $peak_div;                 #Peak division threshold (What precent of max peak?) 

my %final;                    #final information 

my @peaks;                    #Peak numbers
my @peaks_val;                #Peak numbers values
my @eigs;                     #Eigenvalue
my @eigs_raw;                 #Eigenvalue RAW uncut
my @intKS_file;               #intKS File (x - states, y - Raw Projections);
my @intKS_raw_file;           #intKS RAW File (x - states, y - Raw Projections);

my $output = 'peaks.dat';     #main peak output
my $eig_file = 'eig.dat';     #Eigenvalues file
my $gnu = 'temp.gnu';         #gnuplot template file
my $temp = 'output.temp';     #output template file

my ($xr1, $xr2) = (-28, -6);  #for the gnuplot at the end
#--------------------------------------------------------

open my $fh_raw, '<', $file_raw or die "ERROR: Cannot open $file_raw $!";
@intKS_raw_file = <$fh_raw>;

open my $fh, '<', $file or die "ERROR: Cannot open $file $!";
@intKS_file = <$fh>;

close ($fh_raw);
close ($fh);

#first loop through the file find the intKS_file find highest peak
for (@intKS_raw_file){

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
   for (@intKS_raw_file){

      
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
   print $gnu_out "set term x11 title 'Peak division Threshold: $peak_div\n plot \'$file_raw\' u 1:2 w lp, \'$temp\' u 1:2 w p lw 3\n pause -1";
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

         if    ( $ans =~ /y/i)   { last PROMPT; }
         elsif ( $ans =~ /n/i)   { redo FIND_PEAKS;}
         else                    { next PROMPT;}
   }
}

#delete templates
unlink $temp, $gnu;

#Open up the Eigenvalue file
open my $fe, '<', $eig_file or die "ERROR: Cannot Open $eig_file $!";
chomp(@eigs_raw = <$fe>);
close $fe;

#Find the eigenvalues
push @eigs, $eigs_raw[$_ - 1] for (@peaks);

#Find Peak values
for my $val (@eigs){

   #closest value to the eignevaluse, x and y axes
   my ($cur_diff, $cur_peak);

   #loop through the intKS file
   for my $line (@intKS_file){

      my $diff = abs((split ' ', $line)[0] - $val);

      if ( not defined($cur_diff) or ($diff < $cur_diff)  ){
         $cur_peak = (split ' ', $line)[1];
         $cur_diff = $diff;
      }
   }

   #Update the peaks value
   push @peaks_val, $cur_peak;
}


      
#Print the main peaks file
open my $main_out, '>', $output or die "ERROR: Cannot Open $output $!";
select $main_out;

my $ncount;
for my $n (0 .. $#peaks){
   printf " %6d %8d %14.4f %14.7f \n", $n+1, $peaks[$n], $eigs[$n], $peaks_val[$n];
}

#Open gnuplot
my $gnu_temp = 'intKS.gnu';
open my $fh_gnu, '>', $gnu_temp;
select $fh_gnu;

print<<EOF;
set terminal jpeg
set output 'intKS.jpeg'
unset key
set xr[$xr1:$xr2]
set xlabel(\"Electron Binding [eV]\")
set ylabel(\"Photoemission Spectrum [arb. u.]\")
plot './$file' u 1:2 w l, './$output' u 3:(\$4):1 w labels font \"bold,10\"
EOF


system "gnuplot -e \"load \'$gnu_temp\'\"";

select STDIN;
close $fh_gnu;
close $main_out;

