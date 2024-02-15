#include "input.hpp"

void Input::defaultIni(){
    // General parameters
    time_limit = 7200;
    nb_threads = 1;
    seed = 0;
    verbose = 1;
    precision = 0.000001;

    // Instance parameters
    nb_capa = 4;   
    percentage_budget = 0.5;   
    edge_cost_mult = 10;                      

    // Parameters for L-Shaped
    nb_max_it_hot_start = 100;
    papadakos = true;
    use_callback = true;
    user_and_Lazy = false;
    validinequalityEV = true;
    trustregion = false;

    // Parameters for SAA
    nbproblemsSAA = 10;
    nbscenariosSAA = 100;
    nbvalidatescenariosSAA = 1000;
    use_valid_inequality_SAA = true;
    correlated = true;
}

Input::Input(int argc, char* argv[]){
    defaultIni();
    int mandatory = 0;
    for (int i = 1; i < argc; i += 2){
        if(std::string(argv[i]) == "-paramfile"){ // Mandatory parameters
            parameters_file = argv[i+1];
            mandatory += 1;
        }else if(std::string(argv[i]) == "-networkfile"){
            network_file = argv[i+1];
            mandatory += 1;
        }
        else if(std::string(argv[i]) == "-solver"){
            int_solver = std::stoi(argv[i+1]);
            solver = Solver(int_solver);
            mandatory += 1;
        }else if(std::string(argv[i]) == "-distribution"){
            int_distribution = std::stoi(argv[i+1]);
            distribution = Distribution(int_distribution);
            mandatory += 1;
        }
        else if(std::string(argv[i]) == "-timelimit") // Optional : general parameters
            time_limit = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-nb_capa_levels")
            nb_capa = std::stoi(argv[i+1]);  
        else if(std::string(argv[i]) == "-percentage_budget")
            percentage_budget = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-edge_cost_mult")
            edge_cost_mult = std::stod(argv[i+1]);      
        else if(std::string(argv[i]) == "-nbthreads") 
            nb_threads = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-seed") 
            seed = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-verbose")
            verbose = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-precision")
            precision = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-ithotstart") // Optional : L-Shaped solver
            nb_max_it_hot_start = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-papadakos")
            papadakos = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-callback")
            use_callback = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-typecut")
            user_and_Lazy = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-validinequalityEV")
            validinequalityEV = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-trustregion")
            trustregion = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-nbproblemsSAA") // Optional :: SAA problem
            nbproblemsSAA = std::stoi(argv[i+1]); 
        else if(std::string(argv[i]) == "-nbscenariosSAA")
            nbscenariosSAA = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-nbvalidatescenariosSAA")
            nbvalidatescenariosSAA = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-validinequalitySAA")
            use_valid_inequality_SAA = std::stoi(argv[i+1]);
        else{
            std::cerr << "ERROR: Argument '" << argv[i+1] << "' not defined." << std::endl;
            throw std::runtime_error(std::string("Incorrect line of command"));
        }
    }
    defineSolutionFile();
    testParameters();
    if(mandatory != 4){
        std::cerr << "ERROR: Not all mandatory arguments were defined." << std::endl;
        std::cerr << "You need to define: -paramfile -networkfile -solver -distribution";
        throw std::runtime_error(std::string("Incorrect line of command"));
    }
}

void Input::defineSolutionFile(){
    std::string netfile = network_file;
    std::string paramfile = parameters_file;

    while(netfile.find("/") != std::string::npos)
        netfile = netfile.substr(netfile.find("/")+1);
    while(paramfile.find("/") != std::string::npos)
        paramfile = paramfile.substr(paramfile.find("/")+1);
    
    solution_file = "../results/";
    if(solver == Input::Solver::SAA) 
        solution_file += "saa/";
    if(solver == Input::Solver::Linearization || solver == Input::Solver::Linearization_LShaped)
        solution_file += "linearization/";
    if(solver == Input::Solver::Expected ) 
        solution_file += "expected/";
    if(distribution == Input::Distribution::DiscreteChoice)
        solution_file += "discreteselection/";
    if(distribution == Input::Distribution::StdNormalization)
        solution_file += "stdnormalization/";
    if(distribution == Input::Distribution::Normal)
        solution_file += "normal/";
    if(distribution == Input::Distribution::Binomial)
        solution_file += "binomial/";

    
    int cost = edge_cost_mult;
    int perc = percentage_budget*100;
    int seed_aux = seed;

    solution_file += netfile.substr(0, netfile.size()-4) + "_" + paramfile.substr(0, paramfile.size()-4) 
                  + "_budget" + std::to_string(perc) + "_capa" + std::to_string(nb_capa) + "_edgemult" + std::to_string(cost);
    
    if((distribution == Input::Distribution::Binomial || distribution == Input::Distribution::StdNormalization) && (solver == Input::Solver::SAA) ){
        solution_file += "_vi" + std::to_string(use_valid_inequality_SAA);
    }
    if(solver == Input::Solver::Linearization_LShaped){
        solution_file += "_call" + std::to_string(use_callback) + "_ppd" + std::to_string(papadakos) + "_ev" + std::to_string(validinequalityEV);
    }
    if(solver == Input::Solver::SAA){
        solution_file += "_" + std::to_string(nbproblemsSAA);
        solution_file += "_" + std::to_string(nbscenariosSAA);
        solution_file += "_" + std::to_string(nbvalidatescenariosSAA);
    }
    solution_file += "_seed" + std::to_string(seed_aux);
    if(solver == Input::Solver::Linearization || solver == Input::Solver::Linearization_LShaped){
        if(time_limit >= 86400) solution_file += "_long";
    }
    solution_file += ".csv";
    solFile.open(solution_file, std::ios_base::out);
}

void Input::testParameters(){
    if(int_distribution >= 4){
        throw std::runtime_error(std::string("Wrong distribution value."));
        exit(0);
    }
    if(int_solver >= 5){
        throw std::runtime_error(std::string("Wrong solver value."));
        exit(0);
    }
    if(distribution != Input::Distribution::DiscreteChoice && 
            (solver == Input::Solver::Linearization ||solver == Input::Solver::Linearization_LShaped) ){
        throw std::runtime_error(std::string("Linearization problem only defined for discrete choice simple problem."));
        exit(0);
    }
    if(solver == Input::Solver::Expected){
        nbproblemsSAA = 0;
        nbscenariosSAA = 0;
    }
    if(distribution == Input::Distribution::Normal){
        nb_capa = 0;
    }
    if(solver == Input::Solver::Linearization_LShaped){
        if(trustregion == true){
            std::cerr << "WARNING: Trust region not tested properly. Setting parameter to FALSE.\n";
            trustregion = false;
        }
    }
    if(percentage_budget <= -0.0001 || percentage_budget >= 0.9999) {
        std::cerr << "WARNING: Percentage budget must be defined in the interval [0.0,1.0]. It will be set to 0.3 as default.\n";
        percentage_budget = 0.5;
    }
}

void Input::display(){
    std::cout << "------------ PARAMETERS --------------" << std::endl;
    std::cout << "INPUT     FILE: '" << input_file << "'" << std::endl;
	std::cout << "NETWORK   FILE: '" << network_file << "'" << std::endl;
	std::cout << "PARAMETER FILE: '" << parameters_file << "'" << std::endl;
    std::cout << "SOLUTION  FILE: '" << solution_file << "'" << std::endl;
    std::cout << "SEED          :  " << seed << std::endl;
    std::cout << "NB CAPACITY   :  " << nb_capa << std::endl;
    std::cout << "PERC BUDGET   :  " << percentage_budget << std::endl;
    std::cout << "MULT EDGE COST:  " << edge_cost_mult << std::endl;
    std::cout << "SOLVER        :  " << solver << std::endl;
    std::cout << "DISTRIBUTION  :  " << distribution << std::endl;
    std::cout << "TIME LIMIT    :  " << time_limit << std::endl;
    std::cout << "NUMBER THREADS:  " << nb_threads << std::endl;
    std::cout << "VERBOSE       :  " << verbose << std::endl;
    std::cout << "PRECISION     :  " << precision << std::endl;
    if(distribution == Input::Distribution::Binomial || 
            distribution == Input::Distribution::StdNormalization){
        std::cout << "VALID INE SAA :  " << use_valid_inequality_SAA << std::endl;
    }
    if(solver == Input::Solver::Linearization_LShaped){
        std::cout << "IT HOT START  :  " << nb_max_it_hot_start << std::endl;
        std::cout << "PAPADAKOS     :  " << papadakos << std::endl;
        std::cout << "USE CALLBACK  :  " << use_callback << std::endl;
        std::cout << "TYPE CUTS     :  " << user_and_Lazy << std::endl;
        std::cout << "VALID INE EV  :  " << validinequalityEV << std::endl;
        std::cout << "TRUST REGION  :  " << trustregion << std::endl;
    }
    if(solver == Input::Solver::SAA){
        std::cout << "NB PROBLEMS   :  " << nbproblemsSAA << std::endl;
        std::cout << "NB SCENARIOS  :  " << nbscenariosSAA << std::endl;
        std::cout << "NB V SCENARIOS:  " << nbvalidatescenariosSAA << std::endl;
    }
    std::cout << "--------------------------------------" << std::endl;
}

void Input::write(std::vector<std::string> output){
    if(!solFile.fail()){
        for(int i = 0; i < (int)(output.size()/2); ++i){
            solFile << distribution << ";";
            solFile << output[i*2] << ";";
            solFile << network_file << ";";
            solFile << parameters_file << ";";
            solFile << percentage_budget << ";";
            solFile << edge_cost_mult << ";";
            solFile << nb_capa << ";";
            solFile << seed << ";";
            solFile << output[i*2+1]; 
            if(!(output[i*2].find(std::string("SAA")) != std::string::npos)){ // 
                solFile << ";;;;;;;;";
            }
            if(distribution == Input::Distribution::StdNormalization ||
                distribution == Input::Distribution::Binomial){
                solFile << use_valid_inequality_SAA << ";";
            }else{
                solFile << ";";
            }
            if(solver == Input::Solver::Linearization_LShaped){
                solFile << nb_max_it_hot_start << ";";
                solFile << papadakos << ";";
                solFile << use_callback << ";";
                solFile << user_and_Lazy << ";";
                solFile << validinequalityEV << ";";
                solFile << trustregion << ";";
            }else{
                solFile << ";;;;;;";
            }
            solFile << std::endl;
        }
    }
}

void Input::writeHead(){
    if(!solFile.fail()){
        solFile << "DISTRIBUTION;SOLVER;NETWORK_FILE;PARAMETER_FILE;PERC_BUDGET;EDGE_COST_MULT;CAPACITY;SEED;";
        solFile << "BUDGET;TOTAL_NB_SCE;NB_SAA_PROB;NB_SAA_SCE;NB_TEST_SAA;";
        solFile << "STATUS;LB;UB;GAP;TIME;NB_BRANCH_NODES;VAR_BRANCH_NODES;";
        solFile << "STAT_LB;STAT_UB;STAT_GAP;VAR_LB;VAR_UB;VAR_GAP;NB_SAA_NOT_SOLVED;";
        solFile << "VALID_INEQUALITY_SAA;IT_HOT_START;PAPADAKOS;USE_CALLBACK;USER_CUTS;";
        solFile << "VALID_INEQUALITY_EV;TRUST_REGION;"; 
        solFile << std::endl;
    }
}

void Input::help(){
    std::cout << std::endl;
	std::cout << "---------" << std::endl;
	std::cout << std::endl;
}

std::ostream& operator<<(std::ostream& lhs, const Input::Solver & solver) {
    switch(solver) {
        case Input::Solver::SAA: {
            lhs << "SAA";
            break;
        }
        case Input::Solver::Linearization: {
            lhs << "Linearization";
            break;
        }
        case Input::Solver::Linearization_LShaped: {
            lhs << "Linearization LShaped";
            break;
        }
        case Input::Solver::Expected: {
            lhs << "EV";
            break;
        }
        default :{
            lhs << "";
            break;
        }
    }
    return lhs;
}

std::ostream& operator<<(std::ostream& lhs, const Input::Distribution & distribution) {
    switch(distribution) {
        case Input::Distribution::DiscreteChoice: {
            lhs << "Discrete Choice";
            break;
        }
        case Input::Distribution::Binomial: {
            lhs << "Binomial";
            break;
        }
        case Input::Distribution::Normal: {
            lhs << "Normal";
            break;
        }
        case Input::Distribution::StdNormalization: {
            lhs << "Standard Normalization";
            break;
        }
        default :{
            lhs << "";
            break;
        }
    }
    return lhs;
}