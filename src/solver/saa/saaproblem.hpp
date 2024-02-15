#ifndef SAAProblem_HPP
#define SAAProblem_HPP

#include <ctime>
#include <cmath>
#include <iostream>
#include "gurobi_c++.h"
#include "../abstractsolver.hpp"
#include "../../data/saainstance.hpp"
#include "saatestproblem.hpp"
#include "../../data/allsceinstance.hpp"

class SAAProblem : public AbstractSolver {
    protected:
        const SAAInstance *saainstance;
        const Instance *instance;
        SAATest *saatest;

        GRBEnv *env = NULL;
        GRBModel **model = NULL;
        GRBVar **x = NULL;
        GRBVar **y = NULL;
        GRBVar **z = NULL;

        GRBVar **prob = NULL;
        GRBVar ***capa = NULL;
        GRBVar ***cum_capa = NULL;
        GRBVar **mean = NULL;
        GRBVar **stddev = NULL;
        GRBVar ***utility = NULL;
        GRBVar **sum_utility = NULL;
        GRBVar ***mccormick1 = NULL;

        double *x_= NULL;;
        double *y_= NULL;;
        double *prob_= NULL;
        double *mean_= NULL;
        double *stddev_= NULL;
        double **utility_= NULL;
        double *sum_utility_= NULL;

        double feasibility_tol;

        int count_x;
        int count_y;

        double variance_lower;
        double statistical_lower;
        double variance_upper;
        double statistical_upper;
        double variance_gap;
        double statistical_gap;
        bool useinequality;
        double var_nb_branch_nodes;
        int nb_saapr_not_sol_opt;
        double start_time;
        bool exact_problem;

        std::vector<double> lower_estimates;
        std::vector<double> count_bnbnodes;

        void create();
        void defineParameters();
        
        int getNbFoundSolutions(int);
        Status solveProblem(int);
        void updateIteration(Status, int);
        void printInfo(Status, int);
        void deleteVariables();
        void updateVariables(int i);
        void solveTestProblem();
        void computeEstimators();

        void writeSolFile(int);
        void generateDistribution(int);

        void generateCommonVariables(int);
        void generateCommonConstraints(int);
        void discreteChoiceDistribution(int);   
        void normalDistribution(int);
        void binomialDistribution(int);
        void stdNormaDistribution(int);
    
    public:
        SAAProblem(const Input & input);
        ~SAAProblem(); 

        void solve();

        std::vector<std::string> write() const {std::vector<std::string> output;
                                                output.push_back("SAA");
                                                output.push_back(saainstance->write()
                                                    + AbstractSolver::writeline()
                                                    + std::to_string(var_nb_branch_nodes) + std::string(";") 
                                                    + std::to_string(statistical_lower) + std::string(";")  
                                                    + std::to_string(statistical_upper) + std::string(";") 
                                                    + std::to_string(statistical_gap) + std::string(";") 
                                                    + std::to_string(variance_lower) + std::string(";")  
                                                    + std::to_string(variance_upper) + std::string(";") 
                                                    + std::to_string(variance_gap) + std::string(";")
                                                    + std::to_string(nb_saapr_not_sol_opt) + std::string(";"));
                                                return output; }                         
};


#endif