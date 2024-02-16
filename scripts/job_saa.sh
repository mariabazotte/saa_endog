#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=15:00:00
#SBATCH --array=1-800
#SBATCH --output=../slurm/saa-%A_%a.out

i=1

if [ $SLURM_ARRAY_TASK_ID -eq $i ]
then
    cd ../src/
fi

solver=0
verbose=0
timelimit=14400
nbproblems=50
nbscenarios=750
nbvalscenarios=50000
param=../data/params.txt

for seed in 0 1000 2000 3000 4000
do 
    for network in ../data/usa/15nodes4facilities.txt ../data/usa/25nodes4facilities.txt ../data/usa/30nodes4facilities.txt ../data/usa/39nodes4facilities.txt ../data/usa/48nodes4facilities.txt ../data/usa/15nodes5facilities.txt ../data/usa/25nodes5facilities.txt ../data/usa/30nodes5facilities.txt ../data/usa/39nodes5facilities.txt ../data/usa/48nodes5facilities.txt
    do
        for budget in 0.5
        do
            # Normal distribution
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -networkfile $network -paramfile $param -distribution 2 -solver $solver -timelimit $timelimit -verbose $verbose -seed $seed -nb_capa_levels 0 -percentage_budget $budget -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidatescenariosSAA $nbvalscenarios
            fi
            (( i = $i +1 ))

            for capa in 2 3 4
            do
                # Discrete selection
                if [ $SLURM_ARRAY_TASK_ID -eq $i ]
                then
                    ./exe -networkfile $network -paramfile $param -distribution 0 -solver $solver -timelimit $timelimit -verbose $verbose -seed $seed -nb_capa_levels $capa -percentage_budget $budget -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidatescenariosSAA $nbvalscenarios
                fi
                (( i = $i +1 ))

                for d in 1 3
                do
                    for valid in 0 1
                    do
                        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
                        then
                            ./exe -networkfile $network -paramfile $param -distribution $d -solver $solver -timelimit $timelimit -verbose $verbose -seed $seed -nb_capa_levels $capa -percentage_budget $budget -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidatescenariosSAA $nbvalscenarios -validinequalitySAA $valid
                        fi
                        (( i = $i +1 ))
                    done
                done
            done
        done
    done
done

sleep 60