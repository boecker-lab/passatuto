#!/bin/bash

# for db in gnps massbankorbi massbankqtof massbank; do
	# for data in Trees Files; do
		# for method in Original; do
			# qsub -l h_vmem=10G -pe threads 2 -q Sandy.q -o /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method.o -e /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method.e /home/i3scke/MSBlast/callMSBlastDecoyTreesHighMemory.sh -base /home/i3scke/MSBlast/Data -c pos -db $db -d $data -method $method -ppm 10 -ae 2
		# done
	# done
# done

# for db in gnps massbank massbankorbi massbankqtof; do
	# for i in $(seq 10); do
		# for data in Trees Files; do
			# for method in RandomPeaks MixedSpectrum Reroot RandomTree; do
				# qsub -l h_vmem=5G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.o -e /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.e /home/i3scke/MSBlast/callMSBlastDecoyTrees.sh -base /home/i3scke/MSBlast/Data -c pos -db $db -d $data -method $method -addFolder _$i -ppm 10 -ae 2 -rt
			# done
		# done
	# done
# done

# for db in gnps massbank massbankorbi massbankqtof; do
	# for i in $(seq 10); do
		# for data in Trees; do
			# for method in ConditionalFast; do
				# qsub -l h_vmem=5G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.o -e /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.e /home/i3scke/MSBlast/callMSBlastDecoyTrees.sh -base /home/i3scke/MSBlast/Data -c pos -db $db -d $data -method $method -addFolder _$i -ppm 10 -ae 2 -rt
			# done
		# done
	# done
# done

# for db in massbank massbankorbi massbankqtof; do
	# for i in $(seq 10); do
		# for data in Files; do
			# for method in ConditionalFast; do
				# qsub -l h_vmem=5G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.o -e /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.e /home/i3scke/MSBlast/callMSBlastDecoyTrees.sh -base /home/i3scke/MSBlast/Data -c pos -db $db -d $data -method $method -addFolder _$i -ppm 10 -ae 2 -rt
			# done
		# done
	# done
# done

# for db in gnps; do
	# for i in $(seq 10); do
		# for data in Files; do
			# for method in ConditionalFast; do
				# qsub -l h_vmem=10G -pe threads 2 -q Sandy.q -o /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.o -e /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.e /home/i3scke/MSBlast/callMSBlastDecoyTreesHighMemory.sh -base /home/i3scke/MSBlast/Data -c pos -db $db -d $data -method $method -addFolder _$i -ppm 10 -ae 2 -rt
			# done
		# done
	# done
# done

for db in gnps; do
	for i in $(seq 10); do
		for data in Files; do
			for method in MixedSpectrum; do
				qsub -l h_vmem=5G -pe threads 6 -q HighMemory.q -o /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.o -e /home/i3scke/MSBlast/Error/DecoyTrees-$db-$data-$method-$i.e /home/i3scke/MSBlast/callMSBlastDecoyTrees.sh -base /home/i3scke/MSBlast/Data -c pos -db $db -d $data -method $method -addFolder _$i -ppm 10 -ae 2 -rt
			done
		done
	done
done
