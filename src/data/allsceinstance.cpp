#include "allsceinstance.hpp"

bool AllSceInstance::isUncertain(int node, int disr_type) const{
    if(input.getDistribution() == Input::Distribution::Binomial || input.getDistribution() == Input::Distribution::StdNormalization)
        return true;
    else return (node_disr_level[node][disr_type]>0);
}

void AllSceInstance::generateScenarios(){
    total_nb_scenarios = 0;

    // Compute scenarios corresponding to disruption types
    for(int disr=0;disr<nb_disr_type;++disr){
        std::vector<int> capacity;
        addScenarios(capacity,0,disr);
    }

    // Compute scenario corresponding to no disruption
    if(UNCERTAIN_NO_DISR==false){
        std::vector<int> capacity_;
        capacity_.resize(nb_facilities,nb_capa_states);
        scen_facility_level_capa.push_back(capacity_);
        scen_disr_type.push_back(-1);
        ++total_nb_scenarios;
    }
    
    // Initialize scen_facility_capacity (capacity of every facility in each scenario)
    scen_facility_capacity.resize(total_nb_scenarios);
    for(int s = 0; s < total_nb_scenarios; ++s){
        scen_facility_capacity[s].resize(nb_facilities);
        for(int f=0;f<nb_facilities;++f){
            int capa_level = scen_facility_level_capa[s][f];
            scen_facility_capacity[s][f] = capacity[capa_level];
        }
    } 
}

void AllSceInstance::addScenarios(std::vector<int> level_capa, int facility, int disr_type){
    int node = facilities[facility];
    if(facility == nb_facilities-1){
        if(isUncertain(node,disr_type)){
            for(int i=0; i<=nb_capa_states; ++i){
                std::vector<int> real_level_capa(level_capa);
                real_level_capa.push_back(i);
                scen_facility_level_capa.push_back(real_level_capa);
                scen_disr_type.push_back(disr_type);
                ++total_nb_scenarios;
            }
        }else{
            std::vector<int> real_level_capa(level_capa);
            real_level_capa.push_back(nb_capa_states);
            scen_facility_level_capa.push_back(real_level_capa);
            scen_disr_type.push_back(disr_type);
            ++total_nb_scenarios;
        }
    }else{
        if(isUncertain(node,disr_type)){
            for(int i=0; i<=nb_capa_states; ++i){
                std::vector<int> real_level_capa(level_capa);
                real_level_capa.push_back(i);
                addScenarios(real_level_capa,facility+1,disr_type);
            }
        }else{
            level_capa.push_back(nb_capa_states);
            addScenarios(level_capa,facility+1,disr_type);
        } 
    }
}

void AllSceInstance::computeScenarioProbability(const double* x_, const double* prob_, const double** utility_){
    scen_probability.resize(total_nb_scenarios,0.0);
    if(input.getDistribution() == Input::Distribution::DiscreteChoice){
        for(int s = 0; s < total_nb_scenarios; ++s){
            scen_probability[s] = getProbDisruption(s);
            for(int f = 0; f < nb_facilities; ++f){
                for(int p = 0; p < nb_prot_level_facility; ++p){
                    if(x_[f*nb_prot_level_facility+p]>= 0.999){
                        scen_probability[s] = scen_probability[s]*getProbFacility(s,f,p);
                    }
                }
            }
        } 
    }
    else if(input.getDistribution() == Input::Distribution::Binomial){
        std::vector<std::vector<std::vector<double>>> binomial;
        binomial.resize(nb_disr_type);
        for(int dt = 0; dt < nb_disr_type; ++dt){
            binomial[dt].resize(nb_facilities);
            for(int f = 0; f < nb_facilities; ++f){
                for(int capa = 0; capa <= nb_capa_states; ++capa){
                    double comb = std::tgamma(nb_capa_states+1)/(std::tgamma(capa+1)*std::tgamma(nb_capa_states-capa+1));
                    double prob = prob_[dt*nb_facilities+f];
                    double final_prob = std::pow(prob,capa)*std::pow(1.0-prob,nb_capa_states-capa);
                    binomial[dt][f].push_back(comb*final_prob); 
                }
            }
        }  
        for(int s = 0; s < total_nb_scenarios; ++s){
            scen_probability[s] = getProbDisruption(s);
            int dt = scen_disr_type[s];
            if(dt != -1){
                for(int f = 0; f < nb_facilities; ++f){
                    int level = scen_facility_level_capa[s][f];
                    scen_probability[s] = scen_probability[s]*binomial[dt][f][level];
                }
            }
        } 
    }
    else if(input.getDistribution() == Input::Distribution::StdNormalization){
        for(int s = 0; s < total_nb_scenarios; ++s){
            scen_probability[s] = getProbDisruption(s);
            int dt = scen_disr_type[s];
            if(dt != -1){
                for(int f = 0; f < nb_facilities; ++f){
                    int level = scen_facility_level_capa[s][f];
                    scen_probability[s] = scen_probability[s]*utility_[dt*nb_facilities+f][level];
                }
            }
        } 
    }
}

void AllSceInstance::test(){
    if(input.getDistribution() == Input::Distribution::Normal){
        throw std::runtime_error(std::string("AllSeInstance not defined for normal distribution."));
        exit(0);
    }
}