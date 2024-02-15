#ifndef LSHAPED__HPP
#define LSHAPED__HPP

#include <ctime>
#include <iostream>
#include "gurobi_c++.h"
#include "../solver/abstractsolver.hpp"
#include "../data/allsceinstance.hpp"
#include "mainproblems/linearizationmainproblem.hpp"


class LShaped {
    protected:
        const Input & input;
        double start_time;

    public:
        LShaped(const Input & input) : input(input) {}
        virtual ~LShaped() {}
        void solve(AbstractMainProblem*);
};

class LinearizationLShaped : public AbstractSolver {
    private:
        const AllSceInstance *instance = NULL;
        AbstractMainProblem* main = NULL;
        LShaped *lshaped = NULL;
        
    protected:
        void defineParameters(){}
        void create();

    public:
        LinearizationLShaped(const Input & input) : 
                            AbstractSolver(input,"LShaped_Linearization"), 
                            instance(new AllSceInstance(input)), 
                            lshaped(new LShaped(input)) { create();}
        ~LinearizationLShaped() { delete main; delete instance; delete lshaped; }

        void solve();

        std::vector<std::string> write() const {
            std::vector<std::string> output;
            output.push_back("Linearization_LShaped");
            output.push_back(instance->write() + AbstractSolver::writeline());
            return output;
        }
};

#endif