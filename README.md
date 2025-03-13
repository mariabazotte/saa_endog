# saa-endog-opt

## Description

This project contains the code for the paper "Solving Two-Stage Programs with Endogenous
Uncertainty via Random Variable Transformation". It contains the implementation of the two-stage Network Design and Facility Protection Problem (NDFPP) with endogenous uncertainty introduced in this paper.

## Building

### Instances

The folder `data/` contains the data files needed for experiments in the paper. The folder `data/generator` includes the file `usa_network_generator.py` to create new network instances as explained in the paper. The instances are generated based on the Basic US Cities dataset available at https://simplemaps.com/data/us-cities (folder `data/generator/data_usa/`). To generate a network file, execute:

```
python usa_network_generator.py -minpop X -nbfacilities X
```

where:
- -minpop 3000                             # Threshold population to select cities from the southeastern region of the United States. 
- -nbfacilities 4                          # Number of facilities of the instance.

The network file is created at the directory `data/usa/`. The directory `data/` also contains the parameters file, which contains:

- nbProtectionLevelsFacilities=4           # Number of protection levels for the facilities  
- facilitiesImpact=0.3                     # Impact of other facilities (for binomial, normal and std. normalization) 
- penaltyUnmetDemand=10                    # Unit penalty for unmet demand 
- objectiveOrConstriantWithInstallCost=0   # 0: Cost in budget constraint, 1: Cost is in the objective 
- noDisrUncertain=0                        # 0: no-disruption event is not uncertain, 1: no-disruption event is uncertain

### Main Code 

Folder `scr/` contains the source code. To run the source code, execute:

```
make 
./exe -parameters
```

Where we have the following parameters:

#### Mandatory parameters: 
- -paramfile ../data/params.txt                    # Path to parameters file
- -networkfile ../data/usa/15nodes4facilities.txt  # Path to network file 
- -solver 0                                             # 0: SAA Applied to Transformed (Exogenous) Program (SAA-TP in the paper), 1: Det. Equiv. Original (Endogenous) Program (DE-OP in the paper) (Referred here as Linearization), 2: Det. Equiv. Original (Endogenous) Program (DE-OP in the paper) solved with L-Shaped (Referred here as Linearization-LShaped), 3: EV problem and computation of the EEV
- -distribution                                        # 0: Discrete Selection (NDFPP-Selection), 1: Binomial (NDFPP-Binomial), 2: Normal (NDFPP-Normal), 3: Std. Normalization (NDFPP-Discrete)

#### Instance parameters: 
- -nb_capa_levels 2        # Number of capacity levels of facilities for discrete distributions (for Normal is 0)
- -percentage_budget 0.5   # Percentage of maximal budget for the budget constraint 
- -edge_cost_mult  10      # Unit cost for opening one unit of edge 

#### Running parameters:
- -timelimit 3600          # Time limit for optimization 
- -nbthreads 1             # Number of threads for gurobi 
- -seed 0                  # Seed (data generation - cost of facility protection) 
- -verbose 0               # 0 -> nothing, 1-> moderated, 2-> intensive (log files)

#### L-shaped parameters (for DE-OP with L-Shaped): 
- -precision 0.0001     # Precision to stop optimization 
- -ithotstart 100       # Max nb of iterations for hot start 
- -papadakos 1          # 0-> no papadakos cuts, 1-> include papadakos cuts 
- -callback 0           # 0-> normal LShaped (multiple trees), 1-> "Benders and cut" 
- -typecut 0            # 0-> cuts only when a integer solution is found, 1-> cuts for fractional solutions also included 
- -validinequalityEV    # 0-> do not use EV valid inequality, 1-> use EV valid inequality 

#### SAA parameters (for SAA-TP): 
- -nbproblemsSAA 50                 # Number of problems for SAA 
- -nbscenariosSAA 750               # Number of scenarios for each SAA problem 
- -nbvalidatescenariosSAA 150000    # Number of scenarios for the validation problem of the SAA 
- -validinequalitySAA               # 0-> do not use valid inequality for SAA problem, 1-> use valid inequality (only for Binomial and Std. Nomalization)

#### Folders main code

The main code contains the following folders:

- -`src/input` # This folder contains the class Input, which stores information on input parameters provided by the user to the solver during code execution.
- -`src/data`  # This folder contains several classes related to instances: Instance holds general information about an instance, such as network configuration. AllSceInstance enumerates scenarios for the original endogenous program, representing realizations of endogenous variables. SAAInstance generates scenarios for the transformed exogenous program, providing sample realizations of exogenous random variables based on a specific transformation (distribution option).
- -`src/solver` # This folder contains implementations of various approaches. Within `src/solver/saa`, you'll find the SAAProblem class, which implements the SAA-TP approach. In `src/solver/baseline`, there are baseline implementations: the ExpectedValueProblem class handles the EV problem of the original endogenous program, computing the EEV. The LinearizationProblem class implements the DE-OP without L-Shaped.
- `src/lshaped` # This folder contains classes for L-Shaped implementation: Lshaped, LinearizionMainProblem, and LinearizationSubproblem, which implement the DE-OP with L-Shaped.

## Results

The results are saved in a repository such as: 

- results
    - compiled 
    - saa 
        - discreteselection
        - binomial 
        - stdnormalization 
        - normal
    - expected 
        - discreteselection
        - binomial 
        - stdnormalization 
        - normal 
    - linearization: 
        - discreteselection

Each test generates a separate `.csv` file in the folder corresponding to the method—saa (SAA-TP), expected (EV problem), or linearization (DE-OP)—and the endogenous distribution, which can be discrete choice/selection, binomial, stdnormalization, discrete, or normal. The folder `results/compiled` contains an Excel file that compiles the experimental results for the paper. This file includes sheets for each specific instance configuration, providing detailed information about all tests. The data for Table 1 is stored in the `tables_mean` sheet, while Tables 2 and 3 are in `tables_ev_mean` sheet. Figure 1 is based on data from `tables_per_seed` sheet. The table in Appendix K also uses `tables_per_seed` sheet, it considers only instances with seed 0. Similarly, the table in Appendix M is found in `ev_tables_per_seed` sheet and considers instances with seed 0. Finally, the table in Appendix N is contained in `ev_tables_mean` sheet.

## Replicating

The folder `scripts/` contains the sbatch files for replicating the results of the paper. Specifically, the files `job_linearization.sh`, `job_saa.sh`, and `job_expected.sh` contain the commands for the paper's results. After executing these sbatch files, the result `.csv` files are written in the folder `results`, according to the configuration explained before.
        
### Unifying paper results

The file `unify_results.py` unifies the results contained in the `.csv` files, creating the Excel file with the compiled experimental results of the paper. After running the sbatch files, execute:

```
python unify_results.py
```

You can modify this file to see the results of other parameters. The unified results are in the folder `results/compiled`.

## Extra covariance test

We also include an extra test that analyses the covariance generated by using the proposed transformation function. The files are contained in the folder `src/tests`. To execute this, update the make file by changing `main_command.o` to `main_updcorrelation.o`. Include the folder `results/correlation`. Execute:

```
make
./exe -parameters
```
Where we have the following parameters:

### Mandatory parameters:
- -paramfile ../data/params.txt                        # Path to parameters file
- -networkfile ../data/usa/15nodes4facilities.txt      # Path to network file 
- -nbsamples 20                                        # Number of samples
- -sizesamples 1000                                    # Number of scenarios for each sample
- -d 0                                                 # Distribution: 0: Discrete Selection (NDFPP-Selection), 1: Binomial (NDFPP-Binomial), 2: Normal (NDFPP-Normal), 3: Std. Normalization (NDFPP-Discrete)
- -nbsolutions 20                                      # Number of first-stage solutions to evaluate.

By running this code, we obtain the matrix of covariance of the SAA estimates for different first-stage solutions using the same sequence of uniform random numbers (Common Random Numbers). Execute the sbatch file `job_covariance_test.sh` (in folder `scripts/`) for examples.
