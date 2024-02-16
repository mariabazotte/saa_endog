#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --time=16:00:00
#SBATCH --array=1-4
#SBATCH --output=../slurm/correlation-%A_%a.out

i=1

if [ $SLURM_ARRAY_TASK_ID -eq $i ]
then
    cd ../src/
fi

for network in ../data/usa/15nodes4facilities.txt 
do
    for param in ../data/params.txt
    do 
        for d in 0 1 2 3
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -networkfile $network -paramfile $param -nbsamples 20 -sizesamples 1000 -d $d -nbsolutions 20
            fi
            (( i = $i +1 ))
        done
    done
done