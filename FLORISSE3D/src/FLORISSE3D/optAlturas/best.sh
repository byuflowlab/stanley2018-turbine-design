#!/bin/bash

# Passing arguments into the bash script from th$

shear_start=75
shear_end=275
spacing_start=5
spacing_end=15

	for ((i=shear_start; i<=shear_end; i+=100))
	do
		for ((j=spacing_start; j<=spacing_end; j+=5))
		do
			shear=$(echo print $i/1000. | python)
      spacing=$(echo print $j/10. | python)
			awk -v shear=$shear -v spacing=$spacing 'BEGIN {min = 1000000} $1 < min {min=$1} $1 == min {id=FILENAME} {AEP=$2} END {print min, AEP, "#", id}' $shear/$spacing/A_passed/COE* >> min_${shear}.txt
		done
	done

exit 0
