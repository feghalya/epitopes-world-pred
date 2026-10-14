#!/usr/bin/env bash

name=merge_peptides_by_country
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`
echo $name

mkdir -p log output/pickles

fullname=${name}
echo '#!/usr/bin/env bash'"
source /home/feghalya/.virtualenvs/ml/bin/activate
time srun -n 1 ./gisaid/step4b.merge_by_country.py
sleep 10
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --account $RAP_ID --chdir $PWD --time 0-1:00:00 --ntasks 1 --cpus-per-task 8 \
	       --mem-per-cpu 12gb --output log/$fullname.$date.log --error log/$fullname.$date.err --job-name $fullname
