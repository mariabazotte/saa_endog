#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=80G
#SBATCH --time=6:00:00
#SBATCH --array=1-500
#SBATCH --output=../slurm/expected-%A_%a.out

i=1

solver=3
verbose=0
timelimit=14400
nbvalscenarios=150000
param=../instances/params.txt

for seed in 0 1000 2000 3000 4000
do
    for network in ../instances/usa/15nodes4facilities.txt ../instances/usa/25nodes4facilities.txt ../instances/usa/30nodes4facilities.txt ../instances/usa/39nodes4facilities.txt ../instances/usa/48nodes4facilities.txt ../instances/usa/15nodes5facilities.txt ../instances/usa/25nodes5facilities.txt ../instances/usa/30nodes5facilities.txt ../instances/usa/39nodes5facilities.txt ../instances/usa/48nodes5facilities.txt
    do
        for budget in 0.5
        do 
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -networkfile $network -paramfile $param -distribution 2 -solver $solver -timelimit $timelimit -verbose $verbose -seed $seed -nb_capa_levels 0 -percentage_budget $budget -nbvalidatescenariosSAA $nbvalscenarios
            fi
            (( i = $i +1 ))

            for capa in 2 3 4
            do
                for d in 0 1 3
                do
                    if [ $SLURM_ARRAY_TASK_ID -eq $i ]
                    then
                        ./exe -networkfile $network -paramfile $param -distribution $d -solver $solver -timelimit $timelimit -verbose $verbose -seed $seed -nb_capa_levels $capa -percentage_budget $budget -nbvalidatescenariosSAA $nbvalscenarios
                    fi
                    (( i = $i +1 ))
                done
            done
        done
    done
done

sleep 60