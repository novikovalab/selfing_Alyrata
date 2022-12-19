#! /bin/bash
#BSUB -J fastsimcoal[1-100]
#BSUB -q normal
#BSUB -M 600 -R "rusage[mem=1000] "
#BSUB -o fsc_%I_%J_output.txt
#BSUB -e fsc_%I_%J_error.txt


echo "STARTING JOB"
acc=${LSB_JOBINDEX}

model='bottle_scale_asy_mig'
wd='/home/ascott/software/fsc26_linux64/revisions/1mil'

cd ${wd}/${model}
mkdir run_${acc}
cp ${model}* run_${acc}/
cd run_${acc}


fsc26 -t ${model}.tpl  -e ${model}.est -m -n 1000000 -N 1000000 -M 0.01 -l 10 -L 40 -q -c 1 -B 1 -x
