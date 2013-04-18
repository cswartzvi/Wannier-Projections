#! /usr/bin/perl 

use strict;
use warnings;

#Open peak file
my $peak_file = shift @ARGV;
my @origin = (.659635E+01,  0.559085E+01,  0.223083E+01);
my $serial_proj = '/global/homes/c/cswartz/Ext_Programs/Wannier-Projections/proj-serial.x';

#Start, stop and interval for the distance in Angstrom
my $start_radius = 0;
my $dr = 0.1;
my $stop_radius = 6.0;

print "Using origin:\n @origin\n";

open my $fh, '<', $peak_file or die "ERROR: No $peak_file found $!";
my @peak_file = <$fh>;

for (@peak_file){

   my $state = (split ' ', $_)[1];

   #output file
   open my $fout, '>', "output.dat".$state;

   print "State: $state \n";
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


   printf $fout "%7.4f  %7.4f  %7.4f \n",  $radius, $proj, $proj-$prev_proj; 
   printf STDOUT "%7.4f  %7.4f  %7.4f \n",  $radius, $proj, $proj-$prev_proj; 

   $prev_proj = $proj;
   }
 close ($fout);
}
