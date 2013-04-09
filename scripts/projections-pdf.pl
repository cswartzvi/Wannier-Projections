#! /usr/bin/perl
#

use warnings;
use strict;

use PDF::API2::Simple;

my $main_dir = '/home/charles/Desktop/Research/MD_Analysis/MD_water-Ions-64/PBE_novdw_330K/OH/Wannier-Projections/config202100';

print "Main Directory: $main_dir\n";

my $pdf = PDF::API2::Simple->new(
   file => 'output.pdf',
   height => 612,
   width => 792,
);

$pdf->add_font('VerdanaBold');
$pdf->add_font('Verdana');
$pdf->add_page();

my $ldos = $main_dir.'/intKS.jpeg';
$pdf->image($ldos, x => 150, y=> 25, height => 500, width => 550);
$pdf->add_page();

my $ncount;
for my $num (@ARGV){
   $ncount++;

   my $text = "Peak #$ncount ($num)";
   $pdf->next_line();
   $pdf->text($text);
   my $picture_3d = $main_dir.'/Simulation_Files/KS_'.$num.'.png';
   $pdf->image($picture_3d, x => 146, y => 250, width => 630, height => 350);

   #my $picture= $main_dir.'/varying-radius/plot-'.$num.'.jpeg';
   my $picture = $main_dir.'/varying-dist/plot-'.$num.'.jpeg';
   $pdf->image($picture, x => 160, y => 20, width => 550, height => 225);
   $pdf->add_page() unless ($num == $ARGV[$#ARGV]);
}

$pdf->save();
