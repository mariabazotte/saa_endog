#ifndef ABSTRACT_MAIN_PROBLEM_HPP
#define ABSTRACT_MAIN_PROBLEM_HPP

#include <ctime>
#include <iostream>
#include <cmath>
#include "gurobi_c++.h"
#include "../../input/input.hpp"
#include "../../data/instance.hpp"
#include "../../data/saainstance.hpp"
#include "../subproblems/abstractsubproblem.hpp"
#include "../../solver/abstractsolver.hpp"

class AbstractSubProblem;

class AbstractMainProblem : public GRBCallback {
    protected:
        const Input & input;
        AbstractSubProblem *sub = NULL;

        GRBEnv *env = NULL;               /* Variables. */
        GRBModel *model = NULL;
        GRBVar *x = NULL;
        GRBVar *y = NULL;
        GRBVar *theta = NULL;
        GRBVar *mean_z = NULL;
        GRBVar *mean_w = NULL;
        GRBLinExpr ub_callback = 0;

        double lb;
        double ub;
        double gap;
        double nb_branch_nodes;
        Status status;
        double intfeastol;
        double lambda = 0.5;

        int count_x;
        int count_y;
        int count_theta;
        int nb_subproblems;

        double sup = 0.0;
        int nb_via_cuts = 0;
        int nb_opt_cuts = 0;
        int nb_ppd_via_cuts = 0;
        int nb_ppd_opt_cuts = 0;
        int nb_iterations = 0;
        int nb_hot_start_iterations = 0;
        double start_time;
        bool not_init_ppd = true;

        bool papadakos;                     /**< false: Do not include pareto optimal cuts. true: Include it. >*/
        bool callback_;                     /**< false: Do not use callback. true: Use it. >*/
        bool user_cut;                      /**< false: If using callback, do not include user cuts. true: Include it. >*/
        bool valid_inequality_EV;           /**< false: Do not use valid cut (EV). true: Use it. */

        double *x_ = NULL;                  /**< Vectors for first stage solutions. >*/
        double *y_ = NULL; 
        double *theta_ = NULL;

        double *x0_ = NULL;                 /**< Vectors for pareto optimal cuts. (using Papadakos)>*/
        double *y0_ = NULL;

        void defineParameters();
        void createLShapedVariables();

        virtual void callback();
        virtual void initializeVectors();
        virtual void createVariables() = 0;
        virtual void createConstraints() = 0;
        virtual void separate(GRBLinExpr &,int,bool);
    
    public:
        AbstractMainProblem(const Input &, int);
        virtual ~AbstractMainProblem();

        double * getX_() {return x_;}
        double * getY_() {return y_;}

        double * getX0_() {return x0_;}
        double * getY0_() {return y0_;}

        double getLB() const { return lb; }
        double getUB() const { return ub; }
        double getGap() const { return gap; }
        double getNbBranchNodes() const { return nb_branch_nodes; }
        double getIntFeasTol() { return intfeastol; } 
        Status getStatus() const { return status; }
        double getSup() const { return sup; }

        void setStatus(Status n) { status = n; }
        void setSup(double n) { sup = n; }

        int getNbViaCuts() const { return nb_via_cuts; }
        int getNbOptCuts() const { return nb_opt_cuts; }
        int getNbPpdViaCuts() const { return nb_ppd_via_cuts; }
        int getNbPpdOptCuts() const { return nb_ppd_opt_cuts; }
        int getNbIt() const { return nb_iterations; }
        int getNbHotStartIt() const { return nb_hot_start_iterations; }
        double getTimeLimit() const { return std::max(0.0, input.getTimeLimit() - (time(NULL) - start_time)); }

        void increaseNbHotStartIt() { ++nb_hot_start_iterations; }
        void increaseNbIt() { ++nb_iterations; }  
        void setStartTime(double time) { start_time = time; }
        void write_file(std::string name){ model->write(name); }
        void updateUB() { ub = std::min(ub,sup); }
        void updateGAP() { if(ub <= 0.000001 && ub >= -0.000001){ gap = (ub - lb); }else{ gap = (ub - lb) / ub; }}  
        void updateTimeLimit(); 

        void create();
        void solve();
        void solveSubProblems(bool);
        void solveCallback();
        void createTrustRegion();
        void updateTrustRegion();
        void removeTrustRegion(); 
        int getNbFoundSolutions(); 

        virtual void initializeCorePoint();
        virtual void updateCorePoint();
        virtual void updateVariables();
        virtual void deleteVariables();
        virtual void setLinearRelaxation();
        virtual void setInteger(); 
        void verifyVariables();  
};

#endif