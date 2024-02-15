#ifndef ExpectedProblem_HPP
#define ExpectedProblem_HPP


#include "gurobi_c++.h"
#include "../../data/instance.hpp"
#include "../abstractsolver.hpp"
#include "../saa/saatestproblem.hpp"
#include "../../data/saainstance.hpp"
#include "../status.hpp"
#include "../../data/allsceinstance.hpp"


class ExpectedValueProblem : public SimpleAbstractSolver {
    protected:
        const SAAInstance *saainstance;
        const Instance *instance;
        SAATest *saatest;

        GRBVar *x = NULL;
        GRBVar *y = NULL;
        GRBVar *z = NULL;
        GRBVar *w = NULL;

        GRBVar *mean = NULL;
        GRBVar *stddev = NULL;
        GRBVar *prob = NULL;
        GRBVar *aux_expected = NULL;
        GRBVar **A = NULL;

        double *x_ = NULL;
        double *y_ = NULL;

        double *prob_ = NULL;
        double *mean_ = NULL;
        double *stddev_ = NULL;
        double **utility_ = NULL;
        double *sum_utility_ = NULL;

        int count_x;
        int count_y;

        GRBEnv *eev_env = NULL;
        GRBModel *eev_model = NULL;
        GRBVar* eev_z = NULL;

        int count_A;

        double saa_eev_lb;
        double saa_eev_ub;
        double saa_eev_gap;
        double saa_eev_time;
        Status saa_eev_status;
        double saa_eev_stat_ub;
        double saa_eev_var_ub;

        double eev_lb;
        double eev_ub;
        double eev_gap;
        double eev_time;
        Status eev_status;

        int eev_total_nb_scenarios;
        int nb_endogenous_scenarios;

        void create();
        void createDistribution();
        void normalDistribution();
        void binomialDistribution();
        void stdNormaDistribution();
        void discreteChoiceDistribution();
        void integerFeasibility();
        void solveSAAEEV();
        void getVariables();
        void solveEEV();
    
    public:
        ExpectedValueProblem(const Input & input): 
                                        SimpleAbstractSolver(input,"Expected"), saainstance(new SAAInstance(input)),instance(saainstance->getInstance()), 
                                        saatest(new SAATest(input,saainstance)),eev_env(new GRBEnv()), eev_model(new GRBModel(*env))  { 
            defineParameters();
            eev_lb = -GRB_INFINITY;
            eev_ub = GRB_INFINITY;
            eev_gap = 1.0;
            eev_time = 0.0; 

            saa_eev_lb = -GRB_INFINITY;
            saa_eev_ub = GRB_INFINITY;
            saa_eev_gap = 1.0;
            saa_eev_time = 0.0; 

            nb_endogenous_scenarios = instance->getNbDisrType()*std::pow(instance->getNbCapaStates(),instance->getNbFacilities());
            create();
        }
        ~ExpectedValueProblem(){ delete[] x_;
                                 delete[] y_;
                                 if(prob) delete[] prob_;
                                 if(mean) delete[] mean_;
                                 if(stddev) delete[] stddev_;
                                 int nb_disr_type = instance->getNbDisrType();
                                 int nb_facilities = instance->getNbFacilities();
                                 if(aux_expected) { for(int j = 0; j < nb_facilities*nb_disr_type; ++j) delete[] utility_[j]; }
                                 if(aux_expected) delete[] utility_;
                                 if(aux_expected) delete[] sum_utility_;
                                 delete[] x;
                                 delete[] y;
                                 delete[] z;
                                 if(w) delete[] w;
                                 if(mean) delete[] mean;
                                 if(stddev) delete[] stddev;
                                 if(prob) delete[] prob;
                                 if(aux_expected) delete[] aux_expected;
                                 if(A) { for(int i = 0; i < count_A; ++i) delete[] A[i];}
                                 if(A) delete[] A;
                                 if(eev_z)delete[] eev_z;
                                 delete eev_model;
                                 delete eev_env;
                                 delete saatest; 
                                 delete saainstance; }
        void solve();

        std::vector<std::string> write() const {
            std::vector<std::string> output;
            output.push_back("EV");
            output.push_back(";;;;;" + SimpleAbstractSolver::writeline());
            output.push_back("EEV_SAA");
            output.push_back(saainstance->write()
                            + std::to_string(saa_eev_status) + std::string(";")
                            + std::to_string(saa_eev_lb) + std::string(";")
                            + std::to_string(saa_eev_ub) + std::string(";")
                            + std::to_string(saa_eev_gap) + std::string(";")
                            + std::to_string(saa_eev_time) + std::string(";0;0;;")
                            + std::to_string(saa_eev_stat_ub) + std::string(";;;") 
                            + std::to_string(saa_eev_var_ub) + std::string(";;;"));
            if(input.getDistribution() != Input::Distribution::Normal && nb_endogenous_scenarios < 20000){
                output.push_back("EEV");
                output.push_back(std::string(";") + std::to_string(eev_total_nb_scenarios) 
                        + std::string(";;;;")
                        + std::to_string(eev_status) + std::string(";")
                        + std::to_string(eev_lb) + std::string(";")
                        + std::to_string(eev_ub) + std::string(";")
                        + std::to_string(eev_gap) + std::string(";")
                        + std::to_string(eev_time) + std::string(";0;"));
            }
            return output;
        }
};

#endif