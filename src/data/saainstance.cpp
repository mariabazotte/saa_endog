#include "saainstance.hpp"

void SAAInstance::defineTSutudent(){
    if(nb_saa_problems == 1) critical_tstudent = 6.31;
    if(nb_saa_problems == 5) critical_tstudent = 2.02;
    if(nb_saa_problems == 10) critical_tstudent = 1.81;
    else if(nb_saa_problems == 11) critical_tstudent = 1.80;
    else if(nb_saa_problems == 12) critical_tstudent = 1.78;
    else if(nb_saa_problems == 13) critical_tstudent = 1.77;
    else if(nb_saa_problems == 14) critical_tstudent = 1.76;
    else if(nb_saa_problems == 15 || nb_saa_problems == 16) critical_tstudent = 1.75;
    else if(nb_saa_problems == 17) critical_tstudent = 1.74;
    else if(nb_saa_problems == 18 || nb_saa_problems == 19) critical_tstudent = 1.73;
    else if(nb_saa_problems >= 20 && nb_saa_problems <= 22) critical_tstudent = 1.72;
    else if(nb_saa_problems >= 23 && nb_saa_problems <= 26) critical_tstudent = 1.71;
    else if(nb_saa_problems >= 27 && nb_saa_problems <= 39) critical_tstudent = 1.70;
    else if(nb_saa_problems >= 40 && nb_saa_problems <= 59) critical_tstudent = 1.68;
    else if(nb_saa_problems >= 60 && nb_saa_problems <= 79) critical_tstudent = 1.67;
    else if(nb_saa_problems >= 80 && nb_saa_problems <= 199) critical_tstudent = 1.66;
}

void SAAInstance::generateScenarios(){
    cumprob_dtype.resize(instance->getNbDisrType()); // Accumulated probability disruption type
    std::vector<double> prob_disr_type = instance->getProbDisrType();
    std::partial_sum(prob_disr_type.begin(), prob_disr_type.end(), cumprob_dtype.begin());

    // Avoid repeated scenarios 
    nb_sce_not_repeated.resize(nb_saa_problems+1);
    sce_weight.resize(nb_saa_problems+1);

    // Information
    sce_disr_type.resize(nb_saa_problems+1);
    normal_scenario.resize(nb_saa_problems+1);

    if(distribution == Input::DiscreteChoice) generateDiscreteChoiceScenarios();
    else if(distribution == Input::Distribution::Binomial) generateBinomialScenarios();
    else if(distribution == Input::Distribution::Normal) generateNormalScenarios();
    else if(distribution == Input::Distribution::StdNormalization) generateStdNormaScenarios();
}

void SAAInstance::generateDiscreteChoiceScenarios(){
    // Distribution for sampling
    std::uniform_real_distribution<> uniform(0.0,1.0);

    // Parameters for model
    scen_cap_by_prot.resize(nb_saa_problems+1);

    // Accumulated probability capacity level
    std::vector<double> cumprob_capa(instance->getNbCapaStates());

    int num_scenarios = nb_scenarios;
    for(int pr = 0; pr < nb_saa_problems+1; ++pr){
        if(pr == nb_saa_problems) num_scenarios = nb_test_scenarios;
        
        scen_cap_by_prot[pr].resize(num_scenarios);
        nb_sce_not_repeated[pr] = 0;
        normal_scenario[pr] = -1;
        bool updated = false;

        for(int s = 0; s < num_scenarios ; ++s){
            int i = nb_sce_not_repeated[pr];
            int dt = simulateDisruptionType(pr);
            if(dt == -1 && updated == false){
                scen_cap_by_prot[pr][i].resize(instance->getNbFacilities());
                for(int f = 0; f < instance->getNbFacilities(); ++f){
                    scen_cap_by_prot[pr][i][f].resize(instance->getNbFacilityProtLevels(),instance->getFullCapacity());
                }
                updated = true;
            }
            if(dt != -1){
                // Sampling for binomial distribution using uniform (0,1)
                scen_cap_by_prot[pr][i].resize(instance->getNbFacilities());
                for(int f = 0; f < instance->getNbFacilities(); ++f){
                    double uniform_nb = (double) uniform(rd_engine);
                    scen_cap_by_prot[pr][i][f].resize(instance->getNbFacilityProtLevels());
                    for(int p = 0; p < instance->getNbFacilityProtLevels(); ++p){
                        std::vector<double> binomial_prob_capa = instance->getBinomialProbCapaVector(dt,f,p);
                        std::partial_sum(binomial_prob_capa.begin(),
                                        binomial_prob_capa.end(),
                                        cumprob_capa.begin());
                        int arg_capa;
                        arg_capa =  std::distance(cumprob_capa.begin(),
                                        std::upper_bound(cumprob_capa.begin(),
                                                        cumprob_capa.end(), 
                                                        uniform_nb));
                        scen_cap_by_prot[pr][i][f][p] = instance->getCapacity(arg_capa);
                    }
                }
                sce_weight[pr].push_back(1);
                ++nb_sce_not_repeated[pr];
            }
        }
    }
    
    // Compute mean capacity of every facility considering each protection level.
    mean_capa.resize(nb_saa_problems);
    for(int i = 0; i < nb_saa_problems; ++i){
        mean_capa[i].resize(instance->getNbFacilities());
        for(int f = 0; f < instance->getNbFacilities(); ++f){
            mean_capa[i][f].resize(instance->getNbFacilityProtLevels(),0.0);
            for(int p = 0; p < instance->getNbFacilityProtLevels(); ++p){
                mean_capa[i][f][p] = 0.0;
                for(int s = 0; s < nb_sce_not_repeated[i]; ++s){
                    mean_capa[i][f][p] += getWeightScenarioSAA(i,s)*scen_cap_by_prot[i][s][f][p];
                }
            }
        }
    }
}

void SAAInstance::generateBinomialScenarios(){
    // Distribution for sampling
    std::uniform_real_distribution<> uniform(0.0,1.0);

    bernoulli_rv.resize(nb_saa_problems+1);
    ord_bernoulli_rv.resize(nb_saa_problems+1);

    int num_scenarios = nb_scenarios;
    for(int pr = 0; pr < nb_saa_problems+1; ++pr){     
        if(pr == nb_saa_problems) num_scenarios = nb_test_scenarios; // Number scenarios for test is greater

        bernoulli_rv[pr].resize(num_scenarios);
        ord_bernoulli_rv[pr].resize(instance->getNbDisrType());
        nb_sce_not_repeated[pr] = 0;
        normal_scenario[pr] = -1;

        for(int dt = 0; dt < instance->getNbDisrType(); ++dt){
            ord_bernoulli_rv[pr][dt].resize(instance->getNbFacilities());
        }

        for(int s = 0; s < num_scenarios ; ++s){
            int i = nb_sce_not_repeated[pr];
            
            int dt = simulateDisruptionType(pr);
            if(dt != -1){
                // Simulation for bernoulli by taking uniform (0,1)
                bernoulli_rv[pr][i].resize(instance->getNbFacilities());
                for(int f = 0; f < instance->getNbFacilities(); ++f){
                    bernoulli_rv[pr][i][f].resize(instance->getNbBinomialCapaStates());
                    for(int c = 0; c < instance->getNbBinomialCapaStates(); ++c) {
                        bernoulli_rv[pr][i][f][c] = (double)uniform(rd_engine);
                        ord_bernoulli_rv[pr][dt][f].push_back(std::make_tuple(i,c));
                    }
                }
                // Update scenario info
                sce_weight[pr].push_back(1);
                ++nb_sce_not_repeated[pr];
            }
        }

        for(int dt = 0; dt < instance->getNbDisrType(); ++dt){
            for(int n = 0; n < instance->getNbFacilities(); ++n){
                std::sort(ord_bernoulli_rv[pr][dt][n].begin(), ord_bernoulli_rv[pr][dt][n].end(), 
                    [this,pr,n](std::tuple<int,int> & a, std::tuple<int,int> & b){
                        return ( bernoulli_rv[pr][std::get<0>(a)][n][std::get<1>(a)] 
                                < bernoulli_rv[pr][std::get<0>(b)][n][std::get<1>(b)] ); 
                    });
            }
        }
    }
}

void SAAInstance::generateNormalScenarios(){
    // Distribution for sampling
    std::normal_distribution<> normal(0.0,1.0);

    normal_rv.resize(nb_saa_problems+1);
    
    int num_scenarios = nb_scenarios;
    for(int pr = 0; pr < nb_saa_problems+1; ++pr){     
        if(pr == nb_saa_problems) num_scenarios = nb_test_scenarios; // Number scenarios for test is greater

        normal_rv[pr].resize(num_scenarios);
        nb_sce_not_repeated[pr] = 0;
        normal_scenario[pr] = -1;

        for(int s = 0; s < num_scenarios ; ++s){
            int i = nb_sce_not_repeated[pr];

            int dt = simulateDisruptionType(pr);
            if(dt != -1){
                // Simulation for normal by taking standard normal 
                normal_rv[pr][i].resize(instance->getNbFacilities());
                for(int f = 0; f < instance->getNbFacilities(); ++f){
                    normal_rv[pr][i][f] = (double)normal(rd_engine);
                }
                // Update scenario info
                sce_weight[pr].push_back(1);
                ++nb_sce_not_repeated[pr];
            }
        }
    }
    
}

void SAAInstance::generateStdNormaScenarios(){
    // Distributio for sampling 
    std::uniform_real_distribution<> uniform(0.0,1.0);

    stdnorma_rv.resize(nb_saa_problems+1);
    ord_stdnorma_rv.resize(nb_saa_problems+1);

    int num_scenarios = nb_scenarios;
    for(int pr = 0; pr < nb_saa_problems+1; ++pr){     
        if(pr == nb_saa_problems) num_scenarios = nb_test_scenarios; // Number scenarios for test is greater

        stdnorma_rv[pr].resize(num_scenarios);
        ord_stdnorma_rv[pr].resize(instance->getNbDisrType());
        nb_sce_not_repeated[pr] = 0;
        normal_scenario[pr] = -1;

        for(int dt = 0; dt < instance->getNbDisrType() ; ++dt){
            ord_stdnorma_rv[pr][dt].resize(instance->getNbFacilities());
        }

        for(int s = 0; s < num_scenarios ; ++s){
            int i = nb_sce_not_repeated[pr];

            int dt = simulateDisruptionType(pr);
            if(dt != -1){
                // Simulation for softmax by taking uniform (0,1)
                stdnorma_rv[pr][i].resize(instance->getNbFacilities());
                for(int f = 0; f < instance->getNbFacilities(); ++f){
                    stdnorma_rv[pr][i][f] = (double)uniform(rd_engine);
                    ord_stdnorma_rv[pr][dt][f].push_back(i);
                }
                // Update scenario info
                sce_weight[pr].push_back(1);
                ++nb_sce_not_repeated[pr];
            }
        }
        for(int dt = 0; dt < instance->getNbDisrType(); ++dt){
            for(int f = 0; f < instance->getNbFacilities(); ++f){
                std::sort(ord_stdnorma_rv[pr][dt][f].begin(),ord_stdnorma_rv[pr][dt][f].end(), 
                    [this,pr,f](int & a, int & b){
                        return (stdnorma_rv[pr][a][f] < stdnorma_rv[pr][b][f]);
                    });
            }
        }
    }
}

int SAAInstance::simulateDisruptionType(int pr){
    // Simulate the disruption type
    auto upper = std::upper_bound(cumprob_dtype.begin(),cumprob_dtype.end(), (double)rand()/RAND_MAX);

    if(upper == cumprob_dtype.end()){ // No disruption (disruption is not uncertain)
        if(instance->getUncertainNoDisr() == true) {
            std::cerr << "Problem scenario generation. Disruptions should be uncertain." << std::endl;
            throw std::runtime_error(std::string("Incorrect scenario generation."));
            exit(0);
        } 
        if(normal_scenario[pr] == -1){
            // Save the scenario representing no disruption
            normal_scenario[pr] = nb_sce_not_repeated[pr]; 
            ++nb_sce_not_repeated[pr];
            sce_weight[pr].push_back(1);
            sce_disr_type[pr].push_back(-1);
        }else{
            ++sce_weight[pr][normal_scenario[pr]];
        }
        return -1;
    }else{
        // Disruption type of scenario
        int dt = std::distance(cumprob_dtype.begin(), upper);
        sce_disr_type[pr].push_back(dt);
        return dt;
    }
}

double ** SAAInstance::getSolutionCapacityVector(int pr, const double* x_, const double* bern_prob, const double* nor_mean, const double* nor_stddev, const double** stdnorma_prob) const {
    int nb_scenarios = nb_sce_not_repeated[pr];
    int nb_facilities = instance->getNbFacilities();
    int nb_capa_levels = instance->getNbCapaStates();
    double ** capacity = new double*[nb_scenarios];

    if(distribution == Input::Distribution::DiscreteChoice){
        for(int s = 0; s < nb_scenarios; ++s){
            capacity[s] = new double[nb_facilities];
            for(int n = 0; n < nb_facilities; ++n){
                for(int p=0;p< instance->getNbFacilityProtLevels(); ++p)
                capacity[s][n] = std::inner_product(scen_cap_by_prot[pr][s][n].begin(), 
                                                    scen_cap_by_prot[pr][s][n].end(),
                                                    x_+(n*instance->getNbFacilityProtLevels()),0.0);
            }
        }
    }else if(distribution == Input::Distribution::Normal){
        for(int s = 0; s < nb_scenarios; ++s){
            int dt = sce_disr_type[pr][s];
            capacity[s] = new double[nb_facilities];
            for(int n = 0; n < nb_facilities; ++n){
                if(dt != -1){
                    capacity[s][n] = (nor_mean[dt*nb_facilities+n] + normal_rv[pr][s][n]*nor_stddev[dt*nb_facilities+n]);
                }else{
                    capacity[s][n] =  instance->getFullCapacity();
                }
            }
        }
    }else if(distribution == Input::Distribution::Binomial){
        for(int s = 0; s < nb_scenarios; ++s){
            int dt = sce_disr_type[pr][s];
            capacity[s] = new double[nb_facilities];
            for(int n = 0; n < nb_facilities; ++n){
                if(dt != -1){
                    capacity[s][n] = 0.0;
                    for(int w = 0; w < instance->getNbBinomialCapaStates(); ++w){
                        if(bernoulli_rv[pr][s][n][w] < bern_prob[dt*nb_facilities+n]){
                            capacity[s][n] += instance->getCapaPerStage();
                        }
                    }
                }else{
                    capacity[s][n] = instance->getFullCapacity();
                }
            }
        }
    }else if(distribution == Input::Distribution::StdNormalization){
        std::vector<std::vector<std::vector<double>>> cum_stdnorma_prob;
        cum_stdnorma_prob.resize(instance->getNbDisrType());
        for(int dt = 0; dt < instance->getNbDisrType(); ++dt){
            cum_stdnorma_prob[dt].resize(nb_facilities);
            for(int n = 0; n < nb_facilities; ++n){
                cum_stdnorma_prob[dt][n].resize(nb_capa_levels);
                std::partial_sum(stdnorma_prob[dt*nb_facilities+n],
                                 stdnorma_prob[dt*nb_facilities+n]+nb_capa_levels,
                                 cum_stdnorma_prob[dt][n].begin());
            }
        }
        for(int s = 0; s < nb_scenarios; ++s){
            int dt = sce_disr_type[pr][s];
            capacity[s] = new double[nb_facilities];
            for(int n = 0; n < nb_facilities; ++n){
                if(dt != -1){
                    int arg_capa = std::distance(cum_stdnorma_prob[dt][n].begin(),std::upper_bound(
                                                cum_stdnorma_prob[dt][n].begin(),cum_stdnorma_prob[dt][n].end(), 
                                                stdnorma_rv[pr][s][n]));
                    capacity[s][n] = instance->getCapacity(arg_capa);
                }else{
                    capacity[s][n] = instance->getFullCapacity();
                }
            }
        }
    }
    return capacity;
}

