#ifndef LINEARIZATION_SUBPROBLEM_HPP
#define LINEARIZATION_SUBPROBLEM_HPP

#include "gurobi_c++.h"
#include "../../input/input.hpp"
#include "../../data/allsceinstance.hpp"
#include "../mainproblems/abstractmainproblem.hpp"
#include "abstractsubproblem.hpp"
#include "../mainproblems/linearizationmainproblem.hpp"

class LinearizationMainProblem;

class LinearizationSubproblem : public AbstractSubProblem {
    private:
        const AllSceInstance *instance;

        GRBVar *z = NULL;
        GRBVar *w = NULL;
        GRBVar **q = NULL;
        GRBVar ***v = NULL;

        GRBConstr* recursive_one = NULL;
        GRBConstr* recursive_two = NULL;
        GRBConstr* mccormick_one = NULL;
        GRBConstr** mccormick_two = NULL;
        GRBConstr** mccormick_three = NULL;
        GRBConstr** mccormick_four = NULL;
        GRBConstr* arc_capacity = NULL;
        GRBConstr* flow_clients = NULL;
        GRBConstr* flow_facilities = NULL;

        double * u_ = NULL;
        double * v_ = NULL;
        double * l_ = NULL;
        double ** dual_mccormick_three = NULL;
        double ** dual_mccormick_four = NULL;
        
        int count_u;
        int count_v;
        int count_l;
        int count_mccormick;
        int count_mccormick_three;
        int count_mccormick_four;

        void createVariables();
        void createConstraints();

    protected:
        void create();
        void initializeVectors();
        void defineParameters();

    public:
        LinearizationSubproblem(const AllSceInstance *, const Input &);
        ~LinearizationSubproblem();

        double getl_(int a) const { return l_[a]; }
        const double* getv_vector() const { return v_; }
        const double* getu_vector() const { return u_; }
        double getDualMccormickThree(int f,int p,int a) const { return dual_mccormick_three[f][p*instance->getNbArcsWithDummy()+a]; }
        double getDualMccormickFour(int f,int p,int a) const { return dual_mccormick_four[f][p*instance->getNbArcsWithDummy()+a]; }

        void solve(AbstractMainProblem* main, double time, int scenario, bool ppd){
            solve((LinearizationMainProblem*) main, time, scenario, ppd);
        }
        void solve(LinearizationMainProblem*, double, int, bool);
};

#endif