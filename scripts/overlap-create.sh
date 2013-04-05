#! /bin/bash

#First arugment is the proc number
proc=$1
shift


for val in $@
do

#Input Overlap File
infile='overlap.in'${val}
outfile='overlap.out'${val}

cat > $infile << EOF 
&wan_projections
   spec_root    =  'Simulation_Files/WAN_',
   spec_ext     =  '.dat'
   proj_state    =  ${val}
   proj_root     =  'Simulation_Files/KS_',
   proj_ext     =  '.dat'
   nbsp        =  256,
   print_xsf   = .TRUE.
   atomfile    =  'atoms.dat'
   grid(1)     =  128,
   grid(2)     =  128,
   grid(3)     =  128,
   alat(1)     =  23.5170,
   alat(2)     =  23.5170,
   alat(3)     =  23.5170,
   output      =  'overlap.dat${val}'
/
EOF

aprun -n $proc /global/homes/c/cswartz/Ext_Programs/Wannier-Projections/src/wanproj.x < $infile > $outfile

done
