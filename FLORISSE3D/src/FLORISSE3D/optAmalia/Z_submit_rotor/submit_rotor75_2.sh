#!/bin/bash

# Passing arguments into the bash script from th$

spacing_start=5
spacing_end=15
N=100


	for ((k=25; k<=N; k+=1))
	do
		for ((j=spacing_start; j<=spacing_end; j+=5))
		do
			python optLayout_Rotor.py 75 $j 12 $k
		done
	done


exit 0
