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
                       cd ${shear}/${spacing}/
                       good="true"
                       mkdir A_passed
                       mkdir Z_failed
                       while [ "$good" = "true" ];
                       do
                               file=$(echo $(awk '$1=="SNOPTC" && $2=="EXIT" &&  $3=="0" && $5=="finished" && $6=="successfully" {print FILENAME}' summary*.out) | awk -F    '[_.]' '{print $5"_"$6}')
                               if [ "$file" = "_" ]
                               then
                                       mv *.txt Z_failed
                                       mv *.out Z_failed
                                       echo $i $j
                                       good="false"
                               else
                                       mv COE_${file}.txt A_passed
                                       mv diameter_${file}.txt A_passed
                                       mv power_diam_${file}.txt A_passed
                                       mv summary_${file}.out A_passed
                                       mv thickness_${file}.txt A_passed
                                       mv XYZ_${file}.txt A_passed
                                       echo "passed"
                               fi
                       done
                       cd ../..
               done
       done

exit 0
