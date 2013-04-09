#! /usr/bin/perl 

use strict;
use warnings;

#Open peak file
my $peak_file = shift @ARGV;
my @origin = (0.729660E+01,  0.531069E+01,  0.220219E+01);
my $serial_proj = '/home/charles/Desktop/Research/Ext_Programs/Wannier-Projection/proj-serial.x';
print "Using origin:\n @origin\n";

open my $fh, '<', $peak_file or die "ERROR: No $peak_file found $!";
my @peak_file = <$fh>;

for my $state (@peak_file){

   print "State: $state \n";
   #Start, stop and interval for the distance in Angstrom
   my $start_radius = 0;
   my $dr = 0.1;
   my $stop_radius = 6.0;

   #Previous Projection
   my $prev_proj = 0.0;

   for (my $i=0;;$i++){

      my $radius = $i*$dr + $start_radius;
      last if ($radius >= $stop_radius);

      my $infile = 'overlap.in';
      open my $fh, '>', $infile or die "ERROR: Cannot open $infile $!";
      select $fh;

print <<EOF;
&projections
spec_root   =  '../Simulation_Files/KS_',
spec_ext    =  '.dat'
spec_state  =  $state
proj_root   =  '../Simulation_Files/KS_',
proj_ext    =  '.dat'
proj_state  =  $state
print_xsf   = .FALSE.
atomfile    =  'atoms.dat'
grid(1)     =  128,
grid(2)     =  128,
grid(3)     =  128,
alat(1)     =  23.5170,
alat(2)     =  23.5170,
alat(3)     =  23.5170,
int_sphere  =  .TRUE.,
rcut        =  $radius,
defect(1)   =  $origin[0],
defect(2)   =  $origin[1],
defect(3)   =  $origin[2],
/
EOF
   select STDOUT;
   close ($fh);

   my @output = `$serial_proj  < $infile`; 
   my $proj = (split ' ', (grep /Final Projection:/, @output)[0] )[2];


   printf "%7.4f  %7.4f  %7.4f \n",  $radius, $proj, $proj-$prev_proj; 

   $prev_proj = $proj;
   }
}
