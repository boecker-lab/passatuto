#!/bin/bash

base="-base /home/i3scke/MSBlast/searchResultsRPP/ -rpp"

# for db1 in massbankorbi massbankqtof massbank; do
	# for db2 in gnps; do
		# if [ $db1 != $db2 ]; then
			# for method in MassBank CosineDistance OberacherWithIntensities OberacherWithIntensitiesAndMasses; do
				# for data1 in Trees Files; do
					# for data2 in Trees Files; do
						# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-Original.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-Original.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method Original -ppm 10 -ae 2
					# done
				# done
			# done
		# fi
	# done
# done

# for db1 in gnps; do
	# for db2 in massbank; do
		# if [ $db1 != $db2 ]; then
			# for method in MassBank CosineDistance OberacherWithIntensities OberacherWithIntensitiesAndMasses; do
				# for data1 in Trees Files; do
					# for data2 in Trees Files; do
						# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-Original.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-Original.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method Original -ppm 10 -ae 2
					# done
				# done
			# done
		# fi
	# done
# done

# for db1 in massbankorbi massbankqtof massbank; do
	# for db2 in gnps; do
		# if [ $db1 != $db2 ]; then
			# for method in TreeAlignment; do
				# for data1 in Trees; do
					# for data2 in Trees; do
						# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-Original.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-Original.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method Original -ppm 10 -ae 2
					# done
				# done
			# done
		# fi
	# done
# done

# for db1 in gnps; do
	# for db2 in massbank; do
		# if [ $db1 != $db2 ]; then
			# for method in TreeAlignment; do
				# for data1 in Trees; do
					# for data2 in Trees; do
						# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-Original.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-Original.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method Original -ppm 10 -ae 2
					# done
				# done
			# done
		# fi
	# done
# done

# for i in $(seq 10); do
	# for db1 in massbankorbi massbankqtof massbank; do
		# for db2 in gnps; do
			# if [ $db1 != $db2 ]; then
				# for method in MassBank CosineDistance OberacherWithIntensities OberacherWithIntensitiesAndMasses; do
					# for methodDecoy in RandomPeaks MixedSpectrum ConditionalFast Reroot RandomTree; do
						# for data1 in Trees Files; do
							# for data2 in Trees Files; do
								# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -addFolder _$i -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method $methodDecoy -addFolder _$i -ppm 10 -ae 2
							# done
						# done
					# done
				# done
			# fi
		# done
	# done
# done

# for i in $(seq 10); do
	# for db1 in gnps; do
		# for db2 in massbank; do
			# if [ $db1 != $db2 ]; then
				# for method in MassBank CosineDistance OberacherWithIntensities OberacherWithIntensitiesAndMasses; do
					# for methodDecoy in RandomPeaks MixedSpectrum ConditionalFast Reroot RandomTree; do
						# for data1 in Trees Files; do
							# for data2 in Trees Files; do
								# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -addFolder _$i -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method $methodDecoy -addFolder _$i -ppm 10 -ae 2
							# done
						# done
					# done
				# done
			# fi
		# done
	# done
# done

# for i in $(seq 10); do
	# for db1 in massbankorbi massbankqtof massbank; do
		# for db2 in gnps; do
			# if [ $db1 != $db2 ]; then
				# for method in TreeAlignment; do
					# for methodDecoy in Reroot RandomTree; do
						# for data1 in Trees; do
							# for data2 in Trees; do
								# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -addFolder _$i -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method $methodDecoy -addFolder _$i -ppm 10 -ae 2
							# done
						# done
					# done
				# done
			# fi
		# done
	# done
# done

# for i in $(seq 10); do
	# for db1 in gnps; do
		# for db2 in massbank; do
			# if [ $db1 != $db2 ]; then
				# for method in TreeAlignment; do
					# for methodDecoy in Reroot RandomTree; do
						# for data1 in Trees; do
							# for data2 in Trees; do
								# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -addFolder _$i -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method $methodDecoy -addFolder _$i -ppm 10 -ae 2
							# done
						# done
					# done
				# done
			# fi
		# done
	# done
# done

# for i in $(seq 10); do
	# for db1 in massbankorbi massbankqtof massbank gnps; do
		# db2=$db1
		# for method in Equality; do
			# for methodDecoy in RandomPeaks MixedSpectrum ConditionalFast Reroot RandomTree; do
				# for data1 in Trees Files; do
					# data2=$data1
					# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -addFolder _$i -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method $methodDecoy -addFolder _$i -ppm 10 -ae 2
				# done
			# done
		# done
	# done
# done



# for db1 in agilent gnps; do
	# for db2 in agilent gnps; do
		# if [ $db1 != $db2 ]; then
			# for method in MassBank CosineDistance OberacherWithIntensities OberacherWithIntensitiesAndMasses; do
				# for data1 in Trees Files; do
					# for data2 in Trees Files; do
						# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-Original.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-Original.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method Original -ppm 10 -ae 2
					# done
				# done
			# done
		# fi
	# done
# done

# for db1 in agilent gnps; do
	# for db2 in agilent gnps; do
		# if [ $db1 != $db2 ]; then
			# for method in TreeAlignment; do
				# for data1 in Trees; do
					# for data2 in Trees; do
						# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-Original.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-Original.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method Original -ppm 10 -ae 2
					# done
				# done
			# done
		# fi
	# done
# done

# for i in $(seq 10); do
	# for db1 in agilent gnps; do
		# for db2 in agilent gnps; do
			# if [ $db1 != $db2 ]; then
				# for method in MassBank CosineDistance OberacherWithIntensities OberacherWithIntensitiesAndMasses; do
					# for methodDecoy in RandomPeaks MixedSpectrum ConditionalFast Reroot RandomTree; do
						# for data1 in Trees Files; do
							# for data2 in Trees Files; do
								# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -addFolder _$i -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method $methodDecoy -addFolder _$i -ppm 10 -ae 2
							# done
						# done
					# done
				# done
			# fi
		# done
	# done
# done

# for i in $(seq 10); do
	# for db1 in agilent gnps; do
		# for db2 in agilent gnps; do
			# if [ $db1 != $db2 ]; then
				# for method in TreeAlignment; do
					# for methodDecoy in Reroot RandomTree; do
						# for data1 in Trees; do
							# for data2 in Trees; do
								# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -addFolder _$i -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method $methodDecoy -addFolder _$i -ppm 10 -ae 2
							# done
						# done
					# done
				# done
			# fi
		# done
	# done
# done

# for i in $(seq 10); do
	# for db1 in agilent; do
		# db2=$db1
		# for method in Equality; do
			# for methodDecoy in RandomPeaks MixedSpectrum ConditionalFast Reroot RandomTree; do
				# for data1 in Trees Files; do
					# data2=$data1
					# qsub -l h_vmem=4G -pe threads 1 -q Sandy.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -addFolder _$i -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method $methodDecoy -addFolder _$i -ppm 10 -ae 2
				# done
			# done
		# done
	# done
# done

for i in 5 6 7 8 9 10; do
	for db1 in agilent gnps; do
		for db2 in agilent gnps; do
			if [ $db1 != $db2 ]; then
				for method in MassBank CosineDistance OberacherWithIntensities OberacherWithIntensitiesAndMasses; do
					for methodDecoy in ConditionalFast; do
						for data1 in Files; do
							for data2 in Files; do
								qsub -l h_vmem=4G -pe threads 2 -q Default.q -o /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.o -e /home/i3scke/MSBlast/Error/SpectralComparison-$method-$db1-$data1-Original-$db2-$data2-$methodDecoy-$i.e /home/i3scke/MSBlast/callMSBlastSpectralComparison.sh $base -method $method -ppm 10 -ae 2 -postprocess Simple -addFolder _$i -ds1 -base /home/i3scke/MSBlast/Data -c pos -db $db1 -d $data1 -method Original -ppm 10 -ae 2 -ds2 -base /home/i3scke/MSBlast/Data -c pos -db $db2 -d $data2 -method $methodDecoy -addFolder _$i -ppm 10 -ae 2
							done
						done
					done
				done
			fi
		done
	done
done