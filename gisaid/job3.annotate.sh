#!/usr/bin/env bash

name=annotate_sequences
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`
echo $name

mkdir -p log output/scored

for prot in YP_009724389.1 YP_009725295.1 YP_009724390.1 YP_009724391.1 YP_009724392.1 YP_009724393.1 YP_009724394.1 YP_009724395.1 YP_009725318.1 YP_009724396.1 YP_009724397.2 YP_009725255.1
do
    fullname=${name}_${prot}
    echo '#!/usr/bin/env bash'"
source /home/feghalya/.virtualenvs/ml/bin/activate
time srun -n 1 python ./gisaid/step3.annotate.py $prot
sleep 10
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --account $BLUES_ID --workdir $PWD --time 0-6:00:00 --ntasks 1 --cpus-per-task 1 \
	       --mem-per-cpu 64gb --output log/$fullname.$date.log --error log/$fullname.$date.err --job-name $fullname
done
