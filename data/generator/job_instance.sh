#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32

module load gurobi/10.0.3
module load python 
source ~/env_gurobi/bin/activate
pip install utm
pip install pyproj

for minpop in 650000 403500 350000 300000 250000
do
    for nbfacilities in 4 5
    do
        python usa_network_generator.py -minpop $minpop -nbfacilities $nbfacilities
    done
done