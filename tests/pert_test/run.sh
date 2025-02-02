#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
airs=../../src

# Create directory...
rm -rf data && mkdir -p data

# Create perturbation file...
$airs/perturbation data/pert.nc \
		   ../data/AIRS.2003.01.12.166.L1B.AIRS_Rad.v5.0.0.0.G07074102637.hdf \
		   ../data/AIRS.2003.01.12.167.L1B.AIRS_Rad.v5.0.0.0.G07074102639.hdf

# Loop over channel sets...
for pert in 4mu 15mu_low 15mu_high ; do

    # Get map data...
    $airs/map_pert - data/pert.nc data/map_$pert.tab PERTNAME $pert
    
    # Estimate noise...
    $airs/noise_pert - data/pert.nc data/noise_$pert.tab PERTNAME $pert
    
    # Get variance...
    $airs/variance - data/var_$pert.tab data/pert.nc PERTNAME $pert NX 60 NY 30

    # Get events...
    $airs/events - data/events_$pert.tab data/pert.nc PERTNAME $pert VARMIN 0.2
    
done

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.nc data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
