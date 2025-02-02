#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
airs=../../src

# Create directory...
rm -rf data && mkdir -p data

# Extract spectrum...
$airs/spec2tab ../data/AIRS.2003.01.12.166.L1B.AIRS_Rad.v5.0.0.0.G07074102637.hdf index 60 44 data/spec.tab

# Extract quality flags...
../../src/spec_qual ../data/AIRS.2003.01.12.166.L1B.AIRS_Rad.v5.0.0.0.G07074102637.hdf 60 44 data/qual.tab

# Extract map...
$airs/map_rad - \
	      ../data/AIRS.2003.01.12.166.L1B.AIRS_Rad.v5.0.0.0.G07074102637.hdf \
	      ../data/AIRS.2003.01.12.167.L1B.AIRS_Rad.v5.0.0.0.G07074102639.hdf \
	      2338.43 data/wave.tab

# Extract orbit data...
$airs/orbit data/orbit.tab ../data/AIRS.2003.01.12.166.L1B.AIRS_Rad.v5.0.0.0.G07074102637.hdf

exit

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.nc data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
