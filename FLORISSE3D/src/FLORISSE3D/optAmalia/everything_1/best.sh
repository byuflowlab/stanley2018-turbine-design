#!/bin/bash

# Passing arguments into the bash script from th$

shear_start=75
shear_end=275
spacing_start=5
spacing_end=15

rm results*
mkdir XYZ

	for ((i=shear_start; i<=shear_end; i+=100))
	do
		for ((j=spacing_start; j<=spacing_end; j+=5))
		do
			shear=$(echo print $i/1000. | python)
      spacing=$(echo print $j/10. | python)

			file=$(echo $(awk -v shear=$shear -v spacing=$spacing 'BEGIN {min = 1000000} $1 < min {min=$1} $1 == min {id=FILENAME} {AEP=$2} END {print id}' $shear/$spacing/COE*) | awk -F    '[_.]' '{print $4"_"$5}')

			echo SPACING $j >> results_${shear}.txt
			echo >> results_${shear}.txt
			echo file $file >> results_${shear}.txt
			echo coe, aep >> results_${shear}.txt
			echo $(echo $(awk '{print $1, $2}' $shear/$spacing/COE_$file.txt)) >> results_${shear}.txt

			echo rotor diameter >> results_${shear}.txt
			echo $(echo $(awk '{print $2}' $shear/$spacing/power_diam_$file.txt)) >> results_${shear}.txt

			echo hub height >> results_${shear}.txt
			echo $(echo $(awk 'NR==2 {print $3}' $shear/$spacing/XYZ_$file.txt)) >> results_${shear}.txt


			echo rated powers >> results_${shear}.txt
			echo $(echo $(awk '{print $1}' $shear/$spacing/power_diam_$file.txt)) >> results_${shear}.txt


			echo tower diameters >> results_${shear}.txt
			echo $(echo $(awk 'NR==2 {print $1,$2,$3}' $shear/$spacing/diameter_$file.txt)) >> results_${shear}.txt


			echo tower thicknesses >> results_${shear}.txt
			echo $(echo $(awk 'NR==2 {print $1,$2,$3}' $shear/$spacing/thickness_$file.txt)) >> results_${shear}.txt

			echo >> results_${shear}.txt
			echo >> results_${shear}.txt

			cp $shear/$spacing/XYZ_$file.txt XYZ/XYZ_$shear-$spacing.txt
		done
	done

exit 0
