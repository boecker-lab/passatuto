#!/bin/bash


base="-base /home/i3scke/MSBlast/searchResultsRunningTimesRPP/ -rpp"

for i in $(seq 10); do

for db1 in agilent massbankorbi; do
	for db2 in gnps; do
		if [ $db1 != $db2 ]; then
			for method in MassBank CosineDistance OberacherWithIntensities OberacherWithIntensitiesAndMasses; do
				for data1 in Trees Files; do
					for data2 in Trees Files; do
						qsub -l h_vmem=40G -pe threads 12 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparisonRunningTimes-$method-$db1-$data1-Original-$db2-$data2-Original-$i.o -e /home/i3scke/MSBlast/Error/SpectralComparisonRunningTimes-$method-$db1-$data1-Original-$db2-$data2-Original-$i.e /home/i3scke/MSBlast/callMSBlastSpectralComparisonHighMemory.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -rt _$i -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method Original -ppm 10 -ae 2
					done
				done
			done
		fi
	done
done

done