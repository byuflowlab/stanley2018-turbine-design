#!/bin/bash

# Passing arguments into the bash script from th$

shear_start=75
shear_end =275
spacing_start=5
spacing_end=15
N=100


	for ((k=1; k<=N; k+=1))
	do
		for ((j=spacing_start; j<=spacing_end; j+=5))
		do
			python optLayout.py 75 $j 3 $k
		done
	done


exit 0
