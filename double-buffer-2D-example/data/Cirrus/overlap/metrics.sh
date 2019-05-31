#!/bin/bash

FILES=*_3

##Get all required data. The output format is as follow
# LatticeSize NumProc SampleNum TotalTime TimeperLoop TimeHaloExchange
function getData(){

    for f in ${FILES}
    do
	#files are in this fashion LatticeSize_NumProc_SampleNum_API
	echo  $f | awk -F'[_]' ' {printf "%s  %s ", $1, $2}'
#extract Total time and Lattice halo exchange time
	awk '/Loop/ {printf "%1.4f ",$4}; /Computation/ {printf "%1.4f ", $4} ;/Communication/ {printf "%1.4f\n",$4}; ' $f
    done


}

### Start here ###

getData |  sort -n -k 2
