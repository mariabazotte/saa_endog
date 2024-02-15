#include "updcorrelation.hpp"

UpdCorrelation::UpdCorrelation(std::string params_file, std::string network_file, int distribution, 
                                int nb_samples, int size_samples, int nb_solutions):  params_file(params_file),network_file(network_file),input(),
                                                                                      distribution(distribution),nb_samples(nb_samples),size_samples(size_samples),
                                                                                      nb_solutions(nb_solutions){
    create();
}

UpdCorrelation::~UpdCorrelation(){
    for(int i = 0; i < nb_solutions; ++i){
        delete[] x_[i];
        delete[] y_[i];
        if(mean_[i]) delete[] mean_[i];
        if(prob_[i]) delete[] prob_[i];
        int count = instance->getNbFacilities()*instance->getNbDisrType();
        if(utility_[i]) { for(int j = 0; j < count; ++j) if(utility_[i][j]) delete[] utility_[i][j]; }
        if(utility_[i]) delete[] utility_[i];
    }
    delete saatest;
    delete saainstance; 
}

void UpdCorrelation::create(){
    initiateInput();
    saainstance = new SAAInstance(input);
    saatest = new SAATest(input,saainstance);
    instance = (Instance *)saainstance->getInstance();
    nb_facilities = instance->getNbFacilities();
    nb_protec = instance->getNbFacilityProtLevels();
    defineSolutionFile();
}

void UpdCorrelation::initiateInput(){
    input.setDistribution(distribution);
    input.setSolver(Input::Solver::SAA);
    input.setSeed(0);
    input.setNetworkFile(network_file);
    input.setParametersFile(params_file);
    input.setNbProblemsSAA(nb_samples);
    input.setNbScenariosSAA(size_samples);
    input.setNbValidateScenariosSAA(size_samples);
    input.setVerbose(0);
}

void UpdCorrelation::solve(){
    generateSolutions();
    solveSAATest();
    computeAndPrintEstimators();
}

void UpdCorrelation::defineSolutionFile(){
    std::string netfile = input.getNetworkFile();
    std::string paramfile = input.getParametersFile();
    int seed_aux = input.getSeed();
    int budget = 100*input.getPercentageBudget();
    int edgecost = input.getEdgeCostMultiplier();
    int capa = input.getCapaStates();

    while(netfile.find("/") != std::string::npos)
        netfile = netfile.substr(netfile.find("/")+1);
    while(paramfile.find("/") != std::string::npos)
        paramfile = paramfile.substr(paramfile.find("/")+1);

    solution_file = "../results/correlation/";
    if(input.getDistribution() == Input::Distribution::Normal) solution_file += "normal_";
    if(input.getDistribution() == Input::Distribution::Binomial) solution_file += "binomial_";
    if(input.getDistribution() == Input::Distribution::StdNormalization) solution_file += "stdnorma_";
    if(input.getDistribution() == Input::Distribution::DiscreteChoice) solution_file += "discretechoice_";
    
    solution_file += netfile.substr(0, netfile.size()-4) + "_" + paramfile.substr(0, paramfile.size()-4);
    solution_file += "_budget" + std::to_string(budget) + "_capa" + std::to_string(capa) + "_edgemult" + std::to_string(edgecost); 
    solution_file += "_" + std::to_string(nb_samples) + "_" + std::to_string(size_samples) + "_" + std::to_string(nb_solutions);
    solution_file += "_seed" + std::to_string(seed_aux) + ".csv";
}

void UpdCorrelation::solveSAATest(){
    std::cout << "Starting to solve SAA Test." << std::endl;
    objective.resize(nb_solutions);
    for(int s = 0; s < nb_samples; ++s){
        for(int i = 0; i < nb_solutions; ++i){
            double ** capacity_ = saainstance->getSolutionCapacityVector(s,x_[i],prob_[i],mean_[i],stddev_[i],const_cast<const double **>(utility_[i]));
            saatest->solve(x_[i],y_[i],capacity_);
            if(saatest->getStatus() == Status::Optimal) {
                objective[i].push_back(saatest->getUB());
            }else{
                std::cerr << "SAA test is not optimal." << std::endl;
                throw std::runtime_error(std::string("SAA test is not optimal."));
                exit(0);
            }

        }
    }
}

void UpdCorrelation::computeAndPrintEstimators(){
    std::cout << "Printing results." << std::endl; 
    std::ofstream file;
    file.open(solution_file, std::ios_base::out);
    if(!file.fail()){
        file << ";";
        for(int i = 0; i < nb_solutions; ++i) file << "Sol. " << i << ";";
        file << "\n";

        for(int i = 0; i < nb_solutions; ++i){
            file << "Sol. " << i << ";";
            for(int j = 0; j < i; ++j) file << ";";

            double mean_i = std::accumulate(objective[i].begin(), objective[i].end(),0.0)/nb_samples;
            double variance_i = 0.0;
            for(int s = 0; s < nb_samples; ++s) variance_i += std::pow((objective[i][s] - mean_i),2);
            variance_i = variance_i/(nb_samples-1); 
            file << variance_i << ";";

            for(int j = i+1; j < nb_solutions; ++j){
                double mean_j = std::accumulate(objective[j].begin(), objective[j].end(),0.0)/nb_samples;
                double covariance_ij = 0.0;
                for(int s = 0; s < nb_samples; ++s) covariance_ij += (objective[i][s] - mean_i)*(objective[j][s] - mean_j);
                covariance_ij = covariance_ij/(nb_samples-1);
                file << covariance_ij << ";";
            }
            file << "\n";
        }
        file.close();
    }else{
        std::cerr << "ERROR: Unable to open solution file for correlation analysis '" << solution_file << "'." << std::endl; 
        throw std::runtime_error(std::string("Incorrect solution file."));
    }
}

void UpdCorrelation::generateSolutions(){
    if((nb_solutions > 0)){
        randomlyGenerateSolutions();
        // Complete solutions for each distribution
        prob_.resize(nb_solutions,NULL);
        mean_.resize(nb_solutions,NULL);
        stddev_.resize(nb_solutions,NULL);
        utility_.resize(nb_solutions,NULL);
        if(input.getDistribution() == Input::Distribution::Normal) normalDistribution();
        else if(input.getDistribution() == Input::Distribution::Binomial) binomialDistribution();
        else if(input.getDistribution() == Input::Distribution::StdNormalization) stdNormalization();
    }else{
        throw std::runtime_error(std::string("Number of solutions should be greater than 0."));
        exit(0);
    }
}

void UpdCorrelation::randomlyGenerateSolutions(){
    int nb_facilities = instance->getNbFacilities();
    int nb_edges = instance->getNbEdges();
    int nb_protec = instance->getNbFacilityProtLevels();

    x_.resize(nb_solutions);
    y_.resize(nb_solutions);

    std::default_random_engine rd_engine(0);

    std::uniform_int_distribution<> node_protection(0,nb_protec-1);
    std::uniform_int_distribution<> edge_open(0,1);

    for(int i = 0; i < nb_solutions; ++i){
        double budget = 0.0;
        x_[i] = new double[nb_facilities*nb_protec];
        y_[i] = new double[nb_edges];
        for(int f = 0; f < nb_facilities; ++f){
            int selected = (int) node_protection(rd_engine);
            for(int p = 0; p < nb_protec; ++p){
                if(p == selected){
                    x_[i][f*nb_protec+p] = 1;
                    budget += instance->getCostFacilityProtection(f,p);
                }
                else x_[i][f*nb_protec+p] = 0;
            }
        }
        for(int e = 0; e < nb_edges; ++e){
            int open = (int) edge_open(rd_engine);
            if(open == 1 && budget + instance->getCostEdge(e) <= instance->getBudget()){
                y_[i][e] = 1;
                budget += instance->getCostEdge(e);
            }else{
                y_[i][e] = 0;
            }
        }
    }
    std::cout << "Finishing randomly generation of solutions." << std::endl;
}

void UpdCorrelation::binomialDistribution(){
    int count_nodes = instance->getNbFacilities();
    int nb_protec = instance->getNbFacilityProtLevels();
    int nb_disr_type = instance->getNbDisrType();

    double rho = 1 + (instance->getNbFacilities()-1)*instance->getFacilitiesImpact();
    for(int i = 0; i < nb_solutions; ++i){
        prob_[i] = new double[nb_disr_type*count_nodes];
        for(int dt = 0; dt < nb_disr_type; ++dt){
            for(int n = 0; n < count_nodes; ++n){
                prob_[i][dt*count_nodes+n] = 0.0;
                for(int nn = 0; nn < count_nodes; ++nn){
                    double weight = (nn == n) ? 1.0 : instance->getFacilitiesImpact();
                    for(int p = 0; p < nb_protec; ++p){
                        prob_[i][dt*count_nodes+n] += (weight*instance->getBernoulliProbCapa(dt,nn,p)*x_[i][nn*nb_protec+p])/rho;
                    }
                }
            }
        }
    }
}

void UpdCorrelation::normalDistribution(){
    int count_nodes = instance->getNbFacilities();
    int nb_protec = instance->getNbFacilityProtLevels();
    int nb_disr_type = instance->getNbDisrType();
    
    double rho = 1 + (instance->getNbFacilities()-1)*instance->getFacilitiesImpact();
    for(int i = 0; i < nb_solutions; ++i){
        mean_[i] = new double[nb_disr_type*count_nodes];
        stddev_[i] = new double[nb_disr_type*count_nodes];
        for(int dt = 0; dt < nb_disr_type; ++dt){
            for(int n = 0; n < count_nodes; ++n){
                mean_[i][dt*count_nodes+n] = 0.0;
                stddev_[i][dt*count_nodes+n] = 0.0;
                for(int nn = 0; nn < count_nodes; ++nn){
                    double weight = (n == nn) ? 1 : instance->getFacilitiesImpact();
                    stddev_[i][dt*count_nodes+n] += (weight*instance->getStdDevNormalDist(dt,nn))/rho;
                    for(int p = 0; p < nb_protec; ++p){
                        mean_[i][dt*count_nodes+n] += (weight*instance->getMeanNormalDist(dt,nn,p)*x_[i][nn*nb_protec+p])/rho;
                    }
                }
            }
        }
    }
}

void UpdCorrelation::stdNormalization(){
    int count_nodes = instance->getNbFacilities();
    int nb_protec = instance->getNbFacilityProtLevels();
    int nb_disr_type = instance->getNbDisrType();
    int nb_capa_levels = instance->getNbCapaStates();

    for(int i = 0; i < nb_solutions; ++i){
        utility_[i] = new double*[nb_disr_type*count_nodes];
        for(int dt = 0; dt < nb_disr_type; ++dt){
            for(int n = 0; n < count_nodes; ++n){
                double sum_utilities = 0.0;
                utility_[i][dt*count_nodes+n] = new double[nb_capa_levels];
                for(int w = 0; w < nb_capa_levels; ++w){
                    utility_[i][dt*count_nodes+n][w] = 0.0;
                    for(int nn = 0; nn < count_nodes; ++nn){
                        double weight = (n == nn) ? 1 : instance->getFacilitiesImpact();
                        for(int p = 0; p < nb_protec; ++p){
                            utility_[i][dt*count_nodes+n][w] += (weight*instance->getUtility(dt,n,p,w)*x_[i][nn*nb_protec+p]);
                        }
                    }
                    sum_utilities += utility_[i][dt*count_nodes+n][w];
                }
                for(int w = 0; w < nb_capa_levels; ++w){
                    utility_[i][dt*count_nodes+n][w] = utility_[i][dt*count_nodes+n][w]/sum_utilities;
                }
            }
        }
    }
}