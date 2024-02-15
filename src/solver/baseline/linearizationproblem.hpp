#ifndef LinearizationProblem_HPP
#define LinearizationProblem_HPP

#include "gurobi_c++.h"
#include "../abstractsolver.hpp"
#include "../../data/allsceinstance.hpp"

class LinearizationProblem : public SimpleAbstractSolver {
    private:
        const AllSceInstance allinstance;
        GRBVar *x = NULL; 
        GRBVar *y = NULL;
        GRBVar *z = NULL;
        GRBVar *w = NULL;
        GRBVar ***q = NULL;
        GRBVar ****v = NULL;

        double ev_lb;
        double ev_ub;
        double ev_gap;
        double ev_time;
        Status ev_status;
        double eev_lb;
        double eev_ub;
        double eev_gap;
        double eev_time;
        Status eev_status;

        void createVariables();
        void createConstraints();       
    
    protected:
        void defineParameters();
        void create();

    public:
        LinearizationProblem(const Input & input);
        ~LinearizationProblem();
        void solve();
        void solveEEV();

        std::vector<std::string> write() const{
            std::vector<std::string> output;
            output.push_back("Linearization");
            output.push_back(allinstance.write() 
                    + SimpleAbstractSolver::writeline());
            output.push_back("EV");
            output.push_back(allinstance.write() 
                    + std::to_string(ev_status) + std::string(";")
                    + std::to_string(ev_lb) + std::string(";")
                    + std::to_string(ev_ub) + std::string(";")
                    + std::to_string(ev_gap) + std::string(";")
                    + std::to_string(ev_time) + std::string(";0;"));
            output.push_back("EEV");
            output.push_back(allinstance.write() 
                    + std::to_string(eev_status) + std::string(";")
                    + std::to_string(eev_lb) + std::string(";")
                    + std::to_string(eev_ub) + std::string(";")
                    + std::to_string(eev_gap) + std::string(";")
                    + std::to_string(eev_time) + std::string(";0;"));                   
            return output;
        }
};

#endif