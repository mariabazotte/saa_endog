#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

class Input {
    public:
        enum Solver {						
            SAA = 0, 	                 /**< Solved with SAA. **/
            Linearization = 1, 	         /**< Solved with linearization. **/
            Linearization_LShaped = 2,   /**< Solved with LShaped and linearization. **/    
            Expected = 3                 /**< Solving the expected value program. **/
	    };
        
        enum Distribution {
            DiscreteChoice = 0,      /**< Discrete choice of binomials, preprocessing possible. **/
            Binomial = 1,            /**< Function for binomial's parameter (probability of success). **/
            Normal = 2,              /**< Normal distribution approximating binomial. Function for parameters (mean and standard deviation) **/         
            StdNormalization = 3     /**< Ditribution defined by standard normalization function. Function for parameters (utilities) **/
        };

    private:
        const std::string input_file; 
        std::string network_file;      /* File with network information. */
        std::string parameters_file;   /* File with parameters information. */
        std::string solution_file;     /* Directory and file name for solutions. */
        std::ofstream solFile;
        std::string detailed_solution;
        std::ofstream detailedSolFile;

        double time_limit;             /* Time limit for optimization. */
        int nb_threads;                /* Number of threads used in gurobi solver. */
        double seed;                   /* Seed for data generation. */
        int verbose;                   /* Log intensity. */
        double precision;              /* Precision for optimality. */
        int nb_capa;                   /* Nb capa states without 0 capacity (W)*/
        double percentage_budget;      /* Percentage to compute maximum budget. */
        double edge_cost_mult;         /* Edge cost multiplier. */
        bool use_valid_inequality_SAA; /* 0: do not use valid inequality SAA; 1: use it*/
        int nb_max_it_hot_start;       /* Number maximum of hot start iterations. */
        bool papadakos;                /* 0: do not include papadakos cuts, 1: include it. */
        bool use_callback;             /* 0: normal benders, 1: benders using callback. */
        bool user_and_Lazy;            /* 0: only lazy constraints, 1 : lazy constraints and user cuts. */
        bool trustregion;              /* 0: do not use trust region inequalities; 1: include trust region inequalities. */
        bool validinequalityEV;        /* 0: do not include valid inequalities, 1: include valid inequalities. */      
        int nbproblemsSAA;             /* Number of problems used in the SAA. */
        int nbscenariosSAA;            /* Number of scenarios for each SAA problem. */
        int nbvalidatescenariosSAA;    /* Number of scenarios for the validation problem. */
        bool correlated;               /* 0: do not use correlated samples, 1: use correlated samples (only for discrete choice)*/           
        
        int int_solver;                /* Chosen solver for problem. */
        Solver solver;  

        int int_distribution;          /* Chosen distribution. */    
        Distribution distribution;                

        void defaultIni();
        void defineSolutionFile();

    public:
        Input() { defaultIni(); }
        Input(int argc, char* argv[]);
        ~Input(){ solFile.close(); 
                if(!detailedSolFile.fail())  detailedSolFile.close(); }

        std::string getNetworkFile() const { return network_file; }
        std::string getParametersFile() const { return parameters_file; }
        std::string getSolutionFile() const { return solution_file; }
        std::ofstream & getSolFile() { return solFile; }

        // General parameters for the model
        double getTimeLimit() const { return time_limit; }
        int getNbThreads() const { return nb_threads; }
        double getSeed() const { return seed; }
        int getVerbose() const { return verbose; }
        double getPrecision() const { return precision; }

        // Parameters for the instance
        int getCapaStates() const {return nb_capa; }
        double getPercentageBudget() const { return percentage_budget; }
        double getEdgeCostMultiplier() const { return edge_cost_mult; }
        
        // General parameters for L-shaped decomposition
        int getHotStartMaxIt() const { return nb_max_it_hot_start; }
        bool getPapadakos() const { return papadakos; }
        bool getUseCallback() const { return use_callback; }
        bool getUserAndLazy() const { return user_and_Lazy; }
        bool getUseValidInequalityEV() const { return validinequalityEV; }
        bool getTrustRegion() const { return trustregion; }
        
        // General parameters for SAA formulations
        int getNbProblemsSAA() const { return nbproblemsSAA; }
        int getNbScenariosSAA() const { return nbscenariosSAA; }     
        int getNbValidateScenariosSAA() const { return nbvalidatescenariosSAA; }
        bool getUseInequality() const { return use_valid_inequality_SAA; }
        bool getCorrelation() const { return correlated; }
        
        // Definition of the formulation
        const Solver & getSolver() const { return solver; }
        const Distribution & getDistribution() const { return distribution; }

        // Files
        void setNetworkFile(std::string network) { network_file = network; }
        void setParametersFile(std::string param) { parameters_file = param; }
        void setSolutionFile(){ defineSolutionFile(); }

        // Set general parameters
        void setTimeLimit(double time) { time_limit = time; }
        void setNbThreads(int i) { nb_threads = i; }
        void setSeed(double s) { seed = s; }
        void setVerbose(int v) { verbose = v; }
        void setPrecision(double p) { precision = p; }

        // Set info for the instance
        void setCapacity(int capa) { nb_capa = capa; }
        void setUseValidInequalitySAA(bool i) { use_valid_inequality_SAA = i; }
        void setPercentageBudgt(double p) { percentage_budget = p; }
        void setEdgeCostMultiplier(double i) { edge_cost_mult = i; }

		void setItHotStart(int i) { nb_max_it_hot_start = i; }
        void setPapadakos(bool i) { papadakos = i; }
        void setCallback(bool i) { use_callback = i; }
		void setUserCuts(bool i) { user_and_Lazy = i; }
		void setValidInequalityEV(bool i) { validinequalityEV = i; }
        void setTrustRegion(bool i) { trustregion = i; }

        void setNbProblemsSAA(int i) { nbproblemsSAA = i; }
        void setNbScenariosSAA(int i) { nbscenariosSAA = i; }     
        void setNbValidateScenariosSAA(int i) { nbvalidatescenariosSAA = i; }
        void setCorrelation(bool i) { correlated = i; }

        void setSolver(int sol) { int_solver = sol; solver = Solver(int_solver); }
        void setDistribution(int d) { int_distribution = d; distribution = Distribution(int_distribution); }
		
        void testParameters();
        void help();
        void display();
        void write(std::vector<std::string>);
        void writeHead();
};

std::ostream& operator<<(std::ostream& lhs, const Input::Solver & solver);
std::ostream& operator<<(std::ostream& lhs, const Input::Distribution & distribution);


#endif