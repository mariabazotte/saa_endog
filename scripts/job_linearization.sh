#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --array=1-150
#SBATCH --output=../slurm/linearization-%A_%a.out

i=1

if [ $SLURM_ARRAY_TASK_ID -eq $i ]
then
    cd ../src/
fi

timelimit=7200
param=../data/params.txt

for seed in 0 1000 2000 3000 4000
do
    for network in ../data/usa/15nodes4facilities.txt  ../data/usa/25nodes4facilities.txt ../data/usa/30nodes4facilities.txt  ../data/usa/39nodes4facilities.txt ../data/usa/48nodes4facilities.txt ../data/usa/15nodes5facilities.txt  ../data/usa/25nodes5facilities.txt ../data/usa/30nodes5facilities.txt  ../data/usa/39nodes5facilities.txt ../data/usa/48nodes5facilities.txt
    do
        for budget in 0.5 
        do
            for capa in 2 3 4
            do
                if [ $SLURM_ARRAY_TASK_ID -eq $i ]
                then
                    ./exe -networkfile $network -paramfile $param -distribution 0 -solver 2 -timelimit $timelimit -verbose 0 -seed $seed -nb_capa_levels $capa -percentage_budget $budget -callback 0 -papadakos 1 -validinequalityEV 1
                fi
                (( i = $i +1 ))
            done
        done
    done
done

sleep 60