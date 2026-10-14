#!/usr/bin/env bash

name=muscle_alignment
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`
echo $name

ntasks=500

mkdir -p log

echo '#!/usr/bin/env bash'"
module load muscle/3.8.31
source /home/feghalya/.virtualenvs/ml/bin/activate
time srun -n $ntasks python -m mpi4py.futures ./gisaid/step1.parse.py $ntasks
sleep 10
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --account $BLUES_ID --workdir $PWD --time 0-48:00:00 --ntasks $ntasks --cpus-per-task 1 \
	       --mem-per-cpu 8gb --output log/$name.$date.log --error log/$name.$date.err --job-name $name
