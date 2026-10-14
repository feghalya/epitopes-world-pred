#!/usr/bin/env bash

name=netmhc
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`
echo $name

mkdir -p output/netmhc/gisaid log

for pep in 8 9 10 11
do

fullname=${name}.${pep}
echo '#!/usr/bin/env bash'"
module load parallel/20181022 netmhcpan/4.0a
time srun -n 1 parallel -j 32 netMHCpan -p -BA -a {} -l $pep -f data/peptides.$pep.txt -xls -xlsfile output/netmhc/gisaid/{}.$pep.tsv :::: data/alleles.txt > log/$pep.log
sleep 10
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --account $RAP_ID --chdir $PWD --time 0-6:00:00 --ntasks 1 --cpus-per-task 32 --mem-per-cpu 4000mb \
           --output log/$fullname.$date.log --error log/$fullname.$date.err --job-name $fullname

done
