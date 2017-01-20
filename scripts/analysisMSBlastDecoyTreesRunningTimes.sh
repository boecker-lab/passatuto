#!/bin/bash

# for db in gnps; do
	# for i in $(seq 10); do
		# for data in Trees Files; do
			# for method in RandomPeaks MixedSpectrum ConditionalFast Reroot RandomTree; do
				# qsub -l h_vmem=40G -pe threads 12 -q Sandy.q -o /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.o -e /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.e /home/i3scke/MSBlast/callMSBlastDecoyTreesHighMemory.sh -base /home/i3scke/MSBlast/DataRunningTimes -c pos -db $db -d $data -method $method -addFolder _$i -ppm 10 -ae 2 -rt
			# done
		# done
	# done
# done

for db in gnps; do
	for i in $(seq 10); do
		for data in Trees Files; do
			for method in MixedSpectrum Reroot RandomTree; do
				qsub -l h_vmem=40G -pe threads 12 -q Sandy.q -o /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.o -e /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.e /home/i3scke/MSBlast/callMSBlastDecoyTreesHighMemory.sh -base /home/i3scke/MSBlast/DataRunningTimes -c pos -db $db -d $data -method $method -addFolder _$i -ppm 10 -ae 2 -rt
			done
		done
	done
done

