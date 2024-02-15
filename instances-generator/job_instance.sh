#!/cvmfs/soft.computecanada.ca/gentoo/2020/bin/bash
#SBATCH --time=48:00:00
#SBATCH --account=ctb-dionneg1
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32

module load python 
module load gurobi/10.0.3
source ~/env_gurobi/bin/activate
pip install utm
pip install pyproj

python usa_network_generator.py -minpop 300000 -nbfacilities 4