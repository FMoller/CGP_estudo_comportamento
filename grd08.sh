#$ -S /bin/bash

problema=(cm82a)
genmax=(2000000)
tamanho=(100)
semente=(1)
for j in 0
do
	for i in 0
	do
		echo "Solving problem:${problema[j]} with SAM - seed:${semente[i]}"
		#./bin/cgp tables/cm42a.ep seed=1 ncol=100 maxeval=29 mutation=1
		./bin/cgp tables/${problema[j]}.ep seed=${semente[i]} ncol=${tamanho[j]} maxeval=${genmax[j]} mutation=1 resultados/${problema[j]}_${semente[i]} 
	done
done