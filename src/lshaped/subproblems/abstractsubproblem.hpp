#ifndef ABSTRACT_SUBPROBLEM_HPP
#define ABSTRACT_SUBPROBLEM_HPP

#include "gurobi_c++.h"
#include "../../input/input.hpp"
#include "../mainproblems/abstractmainproblem.hpp"
#include "../../solver/abstractsolver.hpp"

class AbstractMainProblem;

class AbstractSubProblem {
    protected:
        const Input & input;

        GRBEnv *env = NULL;
        GRBModel *model = NULL;

        Status status;
        double objective;

        virtual void create() = 0;

        virtual void defineParameters() {
            model->set(GRB_IntParam_Threads, input.getNbThreads());
            model->set(GRB_DoubleParam_TimeLimit, input.getTimeLimit());
            model->set(GRB_IntParam_OutputFlag, 0);
        }

    public:
        AbstractSubProblem(const Input & input) : input(input), env(new GRBEnv()), model(new GRBModel(*env)) {}
        
        virtual ~AbstractSubProblem() {
            delete model;
            delete env; 
        }

        virtual void solve(AbstractMainProblem*, double, int, bool) = 0;

        double getObjective() const { return objective; }
        Status getStatus() const { return status; }
};

#endif