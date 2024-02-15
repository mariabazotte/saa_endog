#ifndef ALLSCEINSTANCE_HPP
#define ALLSCEINSTANCE_HPP

#include "instance.hpp"

class AllSceInstance : public Instance {
    private:
        int total_nb_scenarios;                                       /* Total number of scenarios. */
        std::vector<std::vector<int>> scen_facility_level_capa;       /* Level of capacity facility in each scenario. */
        std::vector<std::vector<int>> scen_facility_capacity;         /* Cpacity of every facility in each scenario; */
        std::vector<int> scen_disr_type;                              /* Disruption type scenario s. */
        std::vector<double> scen_probability;                         /* Probability of scenarios.*/
        void addScenarios(std::vector<int>, int, int);
        
    public:
        AllSceInstance(const Input & input) : Instance(input) { test(); generateScenarios(); }
        ~AllSceInstance() {}

        int getTotalNbScenarios() const { return total_nb_scenarios; }
        double getCapaFacility(int s, int f) const { return scen_facility_capacity[s][f]; } // capacity facility f at scenario s

        bool isUncertain(int,int) const;

        double getProbFacility(int s, int f, int p) const {
            if(scen_disr_type[s] == -1){
                std::cout << "Scenario disruption type -1" << std::endl;
                return 1.0;
            }else{
                int sce_disr = scen_disr_type[s];
                int sce_capa_level = scen_facility_level_capa[s][f];
                return binomial_prob_capa[sce_disr][f][p][sce_capa_level]; 
            }
        }

        double getProbDisruption(int s) const { return scen_disr_type[s] == -1 ? prob_no_disr : prob_disr_type[scen_disr_type[s]] ; }
        double getRhsSimpleCutClient(const double* sol) const { return std::inner_product(client_demand.begin(), client_demand.end(), sol, 0.0000); }
        double getRhsSimpleCutFacility(int s, const double* sol) const { return std::inner_product(scen_facility_capacity[s].begin(), scen_facility_capacity[s].end(), sol, 0.0000); }

        void generateScenarios(); 
        void computeScenarioProbability(const double*, const double*, double const**);
        double getProbabilityScenario(int s) const { return scen_probability[s]; }

        void test();

        std::string write() const { return std::to_string(budget) + std::string(";") + std::to_string(total_nb_scenarios) + std::string(";;;;"); }
};

#endif