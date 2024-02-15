#include "instance.hpp"

Instance::Instance(const Input & input) : input(input), 
                                          network_file(input.getNetworkFile()),
                                          parameters_file(input.getParametersFile()),
                                          seed(input.getSeed()),
                                          nb_capa_states(input.getCapaStates()),
                                          percentage_budget(input.getPercentageBudget()),
                                          max_multi_inst_edge(input.getEdgeCostMultiplier()){
    srand(seed);
    readParametersFile();
    readNetworkFile();
    generateData(); 
    computeProbabilities(); 
}

void Instance::readNetworkFile(){
    std::ifstream file(network_file.c_str());
    if (!file.fail()) {
        file >> nb_nodes;
        file >> nb_edges;
        file >> nb_facilities;
        file >> nb_disr_type; 
        file >> nb_disr_level;

        nb_clients = nb_nodes - nb_facilities;

        std::string name;
        int demand;
        int disr_level;
        int isfacility;
        node_disr_level.resize(nb_nodes);
        for(int n=0; n<nb_nodes; ++n){  
            file >> name;
            file >> demand;
            
            id_cities.insert({name,n});
            name_cities.push_back(name);
            node_demand.push_back(demand);
            for(int i=0; i<nb_disr_type; ++i){
                file >> disr_level;
                int final_disr = ALL_DISR_LEVELS ? disr_level + 1 : disr_level;
                node_disr_level[n].push_back(final_disr);
            }
            
            file >> isfacility;
            if(isfacility){
                facilities.push_back(n);
            }else{
                clients.push_back(n);
                client_demand.push_back(demand);
            }
        } 

        if(nb_facilities != (int)facilities.size() || nb_clients != (int)clients.size()){
            std::string errorMessage = std::string("Error: on file ") + network_file + std::string("\n");
            throw std::runtime_error(errorMessage);
        }

        nb_arcs = 0;
        std::string node1;
        std::string node2;
        double length;
        arcs_leaving_node.resize(nb_nodes);
        arcs_arriving_node.resize(nb_nodes);
        adj_edges.resize(nb_nodes);

        for(int edge=0; edge<nb_edges; ++edge){
            file >> node1;
            file >> node2;
            int n1 = id_cities[node1];
            int n2 = id_cities[node2];

            adj_edges[n1].push_back(edge);
            adj_edges[n2].push_back(edge);

            arcs.push_back({n1,n2});
            edge_of_arc.push_back(edge);
            arcs_leaving_node[n1].push_back(nb_arcs);
            arcs_arriving_node[n2].push_back(nb_arcs);
            ++nb_arcs;

            arcs.push_back({n2,n1});
            edge_of_arc.push_back(edge);
            arcs_leaving_node[n2].push_back(nb_arcs);
            arcs_arriving_node[n1].push_back(nb_arcs);
            ++nb_arcs;

            file >> length;
            edge_length.push_back(length);
        } 

        std::string city_disr;
        double prob_disr;
        for(int i=0; i<nb_disr_type; ++i){
            file >> city_disr;
            file >> prob_disr;
            cities_disr.push_back(city_disr);
            prob_disr_type.push_back(prob_disr);
        }
        prob_no_disr = 1.0 - std::accumulate(prob_disr_type.begin(),prob_disr_type.end(),0.0);

        if(UNCERTAIN_NO_DISR==true){
            ++nb_disr_type;
            prob_disr_type.push_back(prob_no_disr);
            for(int n = 0; n < nb_nodes; ++n )
                node_disr_level[n].push_back(1); // This is used in generating all scenarios
        }

        file.close();
    }else{
        std::string errorMessage = std::string("Error: Could not open file ")+ network_file + std::string("\n");
        throw std::runtime_error(errorMessage);
        exit(0);
    }
}

void Instance::readParametersFile(){
    std::ifstream file(parameters_file.c_str());
    if (file.is_open()) {
        nb_prot_level_facility = std::stoi(getValueFromParameterFile(file,"nbProtectionLevelsFacilities="));
        facilities_impact = std::stod(getValueFromParameterFile(file,"facilitiesImpact="));
        if(facilities_impact <= -0.0001 || facilities_impact >= 0.9999) {
            std::cerr << "WARNING: Facilities impact must be defined in the interval [0.0,1.0]. It will be set to 0.3 as default.\n";
            facilities_impact = 0.3;
        }
        penalty_demand = std::stod(getValueFromParameterFile(file,"penaltyUnmetDemand="));
        const_or_obj_cost = std::stoi(getValueFromParameterFile(file,"objectiveOrConstriantWithInstallCost="));
        if(const_or_obj_cost == true){
            std::cerr << "WARNING: The formulation with cost in the objective has not been tested. It will be set to FALSE (use budget) as default.\n";
            const_or_obj_cost = 0;
        }
        file.close();
    }else {
        std::string errorMessage = std::string("Error: Could not open file ") + parameters_file + std::string("\n");
        throw std::runtime_error(errorMessage);
        exit(0);
    }
}
        
std::string Instance::getValueFromParameterFile(std::ifstream & file, std::string pattern){
    std::string line;
    std::string value = "";
    while(std::getline(file,line)) {
        std::size_t pos = line.find(pattern);
        if (pos != std::string::npos){
            value = line.substr(pos + pattern.size());
            if (value.empty()){
                throw std::runtime_error(std::string("ERROR: Field '" + pattern + "' is empty at parameters file.\n"));
                exit(0);
            }
            return value;
        }
    }
    std::cerr << "WARNING: Did not find '" << pattern << "' inside parameters file." << std::endl; 
    return value;
}

void Instance::generateData(){
    max_flow_arcs = std::accumulate(client_demand.begin(),client_demand.end(),0.0);                                                     
    max_capacity = std::round(max_flow_arcs/(nb_facilities*0.9));                        
    capa_per_state = (double)max_capacity/(double)nb_capa_states;                      
    for(int i = 0; i <= nb_capa_states; ++i){ capacity.push_back(capa_per_state*i); }

    // Cost for protecting facilities
    cost_inst_facility.resize(nb_facilities);
    max_cost_inst_facility.resize(nb_facilities);
    for(int n = 0; n < nb_facilities; ++n){
        max_cost_inst_facility[n] = rand() % 7501 + 7500; // uniform(7500,15000)
        for(int protec = 1; protec <= nb_prot_level_facility; ++protec){
            cost_inst_facility[n].push_back(max_cost_inst_facility[n]*protec/nb_prot_level_facility); 
        }
    }

    // Budget
    double sum_max_budget = 0.0;
    for(int f = 0; f < nb_facilities; ++f){
        sum_max_budget += max_cost_inst_facility[f];
    }
    for(int e = 0; e < nb_edges; ++e){
        sum_max_budget += max_multi_inst_edge*edge_length[e];
    }
    budget = percentage_budget*sum_max_budget;

    // Generate dummy facility 
    double cost_arcs_dummy = (*std::max_element(edge_length.begin(),edge_length.end()))*penalty_demand;
    num_dummy_facility = nb_nodes;
    num_rep_edge_dummy_facility = nb_edges;
    edge_length.push_back(cost_arcs_dummy);
    node_demand.push_back(0);

    arcs_arriving_node.push_back(std::vector<int>()); // we create it, but there is not going to be any arcs arriving at the dummy facility
    arcs_leaving_node.push_back(std::vector<int>()); 

    nb_arcs_dummy = 0;
    for(int i = 0; i < nb_nodes; ++i){
        edge_of_arc.push_back(num_rep_edge_dummy_facility);
        arcs.push_back({num_dummy_facility,i});
        arcs_leaving_node[num_dummy_facility].push_back(nb_arcs+nb_arcs_dummy);
        arcs_arriving_node[i].push_back(nb_arcs+nb_arcs_dummy);
        ++nb_arcs_dummy;
    }
}

void Instance::computeUtilities(){
    // Utilities for standard normalization
    utility.resize(nb_disr_type);
    for(int dt = 0; dt < nb_disr_type; ++dt){
        utility[dt].resize(nb_facilities);
        for(int f = 0; f < nb_facilities; ++f){
            utility[dt][f].resize(nb_prot_level_facility);
            for(int p = 0; p < nb_prot_level_facility; ++p){
                utility[dt][f][p].resize(nb_capa_states+1);
                for(int capa = 0; capa < nb_capa_states+1; ++capa){
                    utility[dt][f][p][capa] = binomial_prob_capa[dt][f][p][capa]*100;
                }
            }
        }
    }
    std::vector<double> sum_maximum;
    sum_maximum.resize(nb_facilities);
    max_utility.resize(nb_disr_type);
    for(int dt = 0; dt < nb_disr_type; ++dt){
        for(int f = 0; f < nb_facilities; ++f){
            double sum_max = 0;
            for(int capa = 0; capa <= nb_capa_states; ++capa){
                double max_val = utility[dt][f][0][capa];
                for(int p = 1; p < nb_prot_level_facility; ++p){
                    if(utility[dt][f][p][capa]>max_val)
                        max_val = utility[dt][f][p][capa];
                }
                sum_max += max_val;
            }
            sum_maximum[f] = sum_max;
        }

        max_utility[dt].resize(nb_facilities);
        for(int f = 0; f < nb_facilities; ++f){
            max_utility[dt][f] = sum_maximum[f];
            for(int f2 = 0; f2 < nb_facilities; ++f2){
                if(f2 != f)
                    max_utility[dt][f] += facilities_impact*sum_maximum[f2];  
            }
        }
    }
}

void Instance::computeFacilitiesMean(){
    mean_capa.resize(nb_facilities);
    for(int n = 0; n < nb_facilities; ++n){
        mean_capa[n].resize(nb_prot_level_facility,0.0);
        if(nb_capa_states > 100){
            for(int p = 0; p < nb_prot_level_facility; ++p){
                mean_capa[n][p] = 0;
                for(int d=0;d<nb_disr_type;++d){
                    mean_capa[n][p] += prob_disr_type[d]*((double)capacity[nb_capa_states])*bernoulli_prob_capa[d][n][p];
                }
                if(UNCERTAIN_NO_DISR==false){
                    mean_capa[n][p] += prob_no_disr*1.0*((double)capacity[nb_capa_states]);
                }
            }
        }else{
            for(int p = 0; p < nb_prot_level_facility; ++p){
                mean_capa[n][p] = 0.0;
                for(int d=0;d<nb_disr_type;++d){
                    for(int capa_state=0;capa_state<=nb_capa_states;++capa_state){
                        mean_capa[n][p] += prob_disr_type[d]*binomial_prob_capa[d][n][p][capa_state]*((double)capacity[capa_state]);
                    }   
                }
                if(UNCERTAIN_NO_DISR==false){
                    mean_capa[n][p] += prob_no_disr*1.0*((double)capacity[nb_capa_states]);
                }
            }
        }
    }
}

void Instance::computeProbabilities(){
    computeFacilitiesProbabilities();
    if(input.getDistribution() == Input::Distribution::DiscreteChoice &&
        (input.getSolver() == Input::Solver::Expected || input.getSolver() == Input::Solver::Linearization || input.getSolver() == Input::Solver::Linearization_LShaped ))
        computeFacilitiesMean();
    if(input.getDistribution() == Input::Distribution::StdNormalization)
        computeUtilities();
    if(input.getDistribution() == Input::Distribution::Normal)
        computeParamsNormalDistr();
}

void Instance::computeFacilitiesProbabilities(){
    bernoulli_prob_capa.resize(nb_disr_type);
    binomial_prob_capa.resize(nb_disr_type);
    
    int nb_disruptions = nb_disr_type;
    if(UNCERTAIN_NO_DISR == true) --nb_disruptions;
    for(int dt = 0; dt < nb_disruptions; ++dt){
        bernoulli_prob_capa[dt].resize(nb_facilities);
        binomial_prob_capa[dt].resize(nb_facilities);
        for(int n = 0; n < nb_facilities; ++n){
            bernoulli_prob_capa[dt][n].resize(nb_prot_level_facility);
            binomial_prob_capa[dt][n].resize(nb_prot_level_facility);
            for(int p = 0; p < nb_prot_level_facility; ++p){
                double disr_level = ALL_DISR_LEVELS ? 
                                    (double)node_disr_level[facilities[n]][dt]/nb_disr_level : 
                                    (double)node_disr_level[facilities[n]][dt]/(nb_disr_level-1); 
                double allocation = (cost_inst_facility[n][p]/max_cost_inst_facility[n])*0.95;   // Multiplier to avoid when the maximum protection is installed, to obtain probability 1 and a single scenario
                double prob = std::pow(allocation,disr_level);
                bernoulli_prob_capa[dt][n][p] = prob;
                for(int capa = 0; capa <= nb_capa_states; ++capa){
                    long double comb = std::tgamma(nb_capa_states+1)/(std::tgamma(capa+1)*std::tgamma(nb_capa_states-capa+1));
                    double final_prob = std::pow(prob,capa)*std::pow(1.0-prob,nb_capa_states-capa);
                    binomial_prob_capa[dt][n][p].push_back(comb*final_prob); 
                }
            }
        }
    }
    if(UNCERTAIN_NO_DISR == true){ // Probabilities are computed slightly different
        bernoulli_prob_capa[nb_disruptions].resize(nb_facilities);
        binomial_prob_capa[nb_disruptions].resize(nb_facilities);
        for(int n = 0; n < nb_facilities; ++n){
            bernoulli_prob_capa[nb_disruptions][n].resize(nb_prot_level_facility);
            binomial_prob_capa[nb_disruptions][n].resize(nb_prot_level_facility);
            for(int p = 0; p < nb_prot_level_facility; ++p){
                double prob = 0.7*(0.95*cost_inst_facility[n][p]/max_cost_inst_facility[n]) + 0.3*(1.0);  
                bernoulli_prob_capa[nb_disruptions][n][p] = prob;
                for(int capa = 0; capa <= nb_capa_states; ++capa){
                    double comb = std::tgamma(nb_capa_states+1)/(std::tgamma(capa+1)*std::tgamma(nb_capa_states-capa+1));
                    double final_prob = std::pow(prob,capa)*std::pow(1.0-prob,nb_capa_states-capa);
                    binomial_prob_capa[nb_disruptions][n][p].push_back(comb*final_prob); 
                }
            }
        }
    }
}

void Instance::computeParamsNormalDistr(){
    mean_normal_dist.resize(nb_disr_type);
    stddev_normal_dist.resize(nb_disr_type);

    double maxP = nb_prot_level_facility;
    double maxL = nb_disr_level;
    double capa = max_capacity;

    int nb_disruptions = nb_disr_type;
    if(UNCERTAIN_NO_DISR == true) --nb_disruptions;
    
    for(int dt = 0; dt < nb_disruptions; ++dt){
        mean_normal_dist[dt].resize(nb_facilities);
        stddev_normal_dist[dt].resize(nb_facilities);
        for(int n = 0; n < nb_facilities; ++n){
            double level = node_disr_level[facilities[n]][dt];
            stddev_normal_dist[dt][n] = (capa*(maxL+level+2.0))/(4.0*2.0*((maxP+1.0)*(maxL+1.0)));
            mean_normal_dist[dt][n].resize(nb_prot_level_facility);
            for(int p = 0; p < nb_prot_level_facility; ++p){
                double interval = 1.0 - ((maxP - (double)p)*((maxL+level+2.0)))/(2.0*(maxP+1.0)*(maxL+1.0));
                mean_normal_dist[dt][n][p] = interval*capa;
            }
        }
    }

    if(UNCERTAIN_NO_DISR == true){
        mean_normal_dist[nb_disr_type-1].resize(nb_facilities);
        stddev_normal_dist[nb_disr_type-1].resize(nb_facilities);
        for(int n = 0; n < nb_facilities; ++n){
            stddev_normal_dist[nb_disr_type-1][n] = (capa*(maxL+2.0))/(4.0*2.0*((maxP+1.0)*(maxL+1.0)));
            mean_normal_dist[nb_disr_type-1][n].resize(nb_prot_level_facility);
            for(int p = 0; p < nb_prot_level_facility; ++p){
                double interval = 1.0 - ((maxP - (double)p)*((maxL+2.0)))/(2.0*(maxP+1.0)*(maxL+1.0));
                mean_normal_dist[nb_disr_type-1][n][p] = interval*capa;
            }
        }
    }
}

std::string Instance::write() const{
    std::string output = std::string("");
    output += std::to_string(nb_nodes) + std::string(";");
    output += std::to_string(nb_edges) + std::string(";");
    output += std::to_string(nb_facilities) + std::string(";");
    output += std::to_string(nb_prot_level_facility) + std::string(";");
    output += std::to_string(nb_capa_states) + std::string(";");
    if(input.getDistribution() == Input::Distribution::Binomial || input.getDistribution() == Input::Distribution::StdNormalization){
        output += std::to_string(facilities_impact) + std::string(";");
    }else{
        output += std::string("-;");
        output += std::string("-;");
    }
    output += std::to_string(penalty_demand) + std::string(";");
    output += std::to_string(const_or_obj_cost) + std::string(";");
    if(const_or_obj_cost){
        output += std::string("-;");
    }else{
        output += std::to_string(budget) + std::string(";");
    }
    return output;
}