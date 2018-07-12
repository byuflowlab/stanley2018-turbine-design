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
                    	  awk -v shear=$shear -v spacing=$spacing 'BEGIN {min = 1000000} $1 < min {min=$1} $1 == min {id=FILENAME} {AEP=$2} END {print min, AEP, "#", id}' $shear/$spacing/COE* >> min_${shear}.txt
                done
       done

       for ((i=shear_start; i<=shear_end; i+=100))
       do
               for ((j=spacing_start; j<=spacing_end; j+=5))
               do
                       shear=$(echo print $i/1000. | python)
                       spacing=$(echo print $j/10. | python)
                       line=$((j / 5))
                       #  echo $line
                       COE=$(echo $(awk -v line="$line" 'NR==line {print $1}' min_$shear.txt))
                       #  echo $COE
                       cd ${shear}/${spacing}/
                       file=$(echo $(awk -v COE="$COE" '$1==COE {print FILENAME}' COE*.txt) | awk -F    '[_.]' '{print $2"_"$3}')
                      #  echo $file
                       p=$(echo $(awk -v file="$file" '$1=="SNOPTC" && $2=="EXIT" &&  $3=="0" && $5=="finished" && $6=="successfully" {print "pass"}' summary_$file.out))
                       echo $p
                       if [ $p == "pass" ]
                       then
                         pass="pass"
                       else
                         pass="fail"
                       fi
                       echo SHEAR: $shear, SPACING: $spacing, FILE: $file: $pass
                       cd ../..
               done
       done

       rm min*

exit 0
