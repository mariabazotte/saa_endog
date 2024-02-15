#ifndef SAAINSTANCE_HPP
#define SAAINSTANCE_HPP

#include "instance.hpp"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <random>

class SAAInstance {
    protected:
        const Instance* instance;                       /* Instance info for scenario generation. */
        Input::Distribution distribution;               /* Distribution. */
        double seed;                                    /* Seed for random data generation. */
        std::default_random_engine rd_engine;           /* Random engine for data generation.*/  
        int nb_saa_problems;                            /* Number of problems for SAA. */
        int nb_scenarios;                               /* Number of scenarios for SAA. */
        int nb_test_scenarios;                          /* Number of scenarios to test the best SAA problem solution. */
        double critical_tstudent;                       /* Critical value for t-student for the nb of problems - according to 95%.*/
        
        std::vector<int> nb_sce_not_repeated;           /* Number of scenarios for each problem (including test scenario) of the SAA without repeated scenarios. */
        std::vector<std::vector<double>> sce_weight;    /* Weight scenario for each problem in the objective for SAA. */ 
        std::vector<int> normal_scenario;               /* Scenario representing no disruption for each problem. */
        std::vector<std::vector<double>> sce_disr_type; /* Disruption type for each scenario. */
        std::vector<double> cumprob_dtype;              /* Accumulated probability of disruption types. */

        std::vector<std::vector<std::vector<std::vector<double>>>> scen_cap_by_prot;              /* scen_cap_by_prot[pr][s][f][p]   : Capacity of facility f with protection level p in scenario s for each problem pr. */
        std::vector<std::vector<std::vector<std::vector<double>>>> bernoulli_rv;                  /* bernoulli_rv[pr][s][f][c]       : Random sample from U(0,1) for problem pr, scenario s, facility f, capacity level c. */
        std::vector<std::vector<std::vector<std::vector<std::tuple<int,int>>>>> ord_bernoulli_rv; /* ord_bernoulli_rv[pr][dt][f][i]  : Ordered list (considering random sample) of indexes tuple<scenario,capacity> for problem pr, disruption dt and facility f.*/
        std::vector<std::vector<std::vector<double>>> stdnorma_rv;                                /* stdnorma_rv[pr][s][f]           : Random sample from U(0,1) for problem pr, scenario s, facility f. */
        std::vector<std::vector<std::vector<std::vector<int>>>> ord_stdnorma_rv;                  /* ord_stdnorma_rv[pr][dt][f][i]   : Ordered list (considering random sample) of indexes vector<scenario> for problem pr, disruption dt and facility f.*/
        std::vector<std::vector<std::vector<double>>> normal_rv;                                  /* normal_rv[pr][s][f]             : Random sample from std normal for problem pr, scenario s, facility f. */
        std::vector<std::vector<std::vector<double>>> mean_capa;                                  /* mean_capa[pr][f][p]             : Mean capacity of facility f considering protection level p in problem pr.*/


        void generateScenarios();                            /* General function to create scenarios. */
        int simulateDisruptionType(int);                     /* Function to simulate disruption type of a scenario. */
        void generateBinomialScenarios();                    /* Scenarios for binomial distribution defined by a function in the formulation.*/
        void generateNormalScenarios();                      /* Scenarios for normal distribution. */
        void generateStdNormaScenarios();                    /* Scenarios for std normalization distribution. */
        void generateDiscreteChoiceScenarios();              /* Scenarios for discrete choice distribution. */
        void defineTSutudent();                              /* Define the critical value of the t-student according to the nb of problems. */
    
    public:
        SAAInstance(const Input & input): instance(new Instance(input)), distribution(input.getDistribution()),seed(0), rd_engine(seed),
                                          nb_saa_problems(input.getNbProblemsSAA()), nb_scenarios(input.getNbScenariosSAA()),
                                          nb_test_scenarios(input.getNbValidateScenariosSAA()){ 
            defineTSutudent();
            generateScenarios();
        }
       
        virtual ~SAAInstance() { delete instance; }

        const Instance * getInstance() const { return instance; }

        int getNbProblems() const { return nb_saa_problems; }
        long long getOriginalNbTestScenarios() const { return nb_test_scenarios; }
        int getNbScenarios(int problem) const { return nb_sce_not_repeated[problem]; }
        int getNoDisruptionScenario(int problem) const { return normal_scenario[problem]; }
        int getNbTestScenarios() const { return nb_sce_not_repeated[nb_saa_problems]; }
        int getTestProblem() const { return nb_saa_problems; }
        double getWeightScenarioSAA(int problem, int sce) const { return sce_weight[problem][sce]/nb_scenarios; }
        double getWeightTestScenarios(int sce) const { return sce_weight[nb_saa_problems][sce]/nb_test_scenarios; }
        double getNbRepeatedTestScenarios(int sce) const { return sce_weight[nb_saa_problems][sce]; }
        double getScenarioDisruption(int problem, int sce) const { return sce_disr_type[problem][sce]; }
        Input::Distribution getDistribution() const { return distribution; }
        double getCriticalTStudent() const { return critical_tstudent; }
        double getCriticalNormal() const { return 1.64; }

        double getCapaFacilityProtection(int problem, int s, int f, int p) const { return scen_cap_by_prot[problem][s][f][p]; }
        double getMeanCapacity(int pr,int f,int p) const { return mean_capa[pr][f][p]; }
        double getNormalDistRV(int pr, int s, int f) const { return normal_rv[pr][s][f]; }
        double getBinomialUniformRV(int pr, int s, int f, int c) const { return bernoulli_rv[pr][s][f][c]; }
        double getStdNormaUniformRV(int pr, int s, int f) const { return stdnorma_rv[pr][s][f]; }
        const std::vector<std::tuple<int,int>> & getOrderedBernoulliUniformRV(int pr, int dt, int f) const { return ord_bernoulli_rv[pr][dt][f]; }
        const std::vector<int> & getOrderedStdNormaUniformRV(int pr, int dt, int f) const { return ord_stdnorma_rv[pr][dt][f]; }
        double ** getSolutionCapacityVector(int pr, const double* x_, const double* bern_prob, const double* nor_mean, const double* nor_stddev, double const** soft_prob) const;

        std::string write() const{
            return std::to_string(instance->getBudget()) + std::string(";") +
                   std::string(";") + std::to_string(nb_saa_problems) + 
                   std::string(";") + std::to_string(nb_scenarios) +
                   std::string(";") + std::to_string(nb_test_scenarios) + std::string(";");
        }
};

#endif