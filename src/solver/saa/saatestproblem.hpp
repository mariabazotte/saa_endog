#ifndef SAA_TEST_HPP
#define SAA_TEST_HPP

#include "gurobi_c++.h"
#include "../../input/input.hpp"
#include "../abstractsolver.hpp"
#include "../../data/saainstance.hpp"
#include "../../lshaped/mainproblems/abstractmainproblem.hpp"
#include <vector>
#include <iostream>

class SAATest : public AbstractSolver {
   protected:
      const SAAInstance *saainstance;
      const Instance *instance;
      int num_scenarios;
      
      GRBEnv *env = NULL;
      GRBModel **model = NULL;
      GRBVar **z = NULL;
      GRBVar **w = NULL;

      AllSceInstance *allsceinstance = NULL;
      GRBEnv *exact_env = NULL;
      GRBModel *exact_model = NULL;
      GRBVar *exact_z = NULL;
      bool exact_problem;
      
      void create();
      void defineParameters();
      void solve();
      virtual void createProblem(int);
      void createExactTestProblem();
      void solveExactTestProblem(double *,double *);

      std::vector<double> cost_scenarios; /* Cost of each scenario in the optimal solution. */

   public:
      SAATest(const Input & input,const SAAInstance* saainstance);
      virtual ~SAATest() { for(int i = 0; i < num_scenarios; ++i){
                              delete[] z[i];
                              delete[] w[i];
                              delete model[i];
                           }
                           delete[] z;
                           delete[] w;
                           delete[] model;
                           delete env;

                           if(allsceinstance) delete allsceinstance;
                           if(exact_z) delete[] exact_z;
                           if(exact_model) delete exact_model;
                           if(exact_env) delete exact_env; 
                        }
      
      void setExactProblem(bool i) { exact_problem = i;}
      void solve(double *,double *,double **);
      void solve(double *,double *,double*,double*,double*,double **);
      double getCostScenario(int s) const { return cost_scenarios[s]; }
      void updateBounds(int s){
         lb += model[s]->get(GRB_DoubleAttr_ObjBound);
         ub += model[s]->get(GRB_DoubleAttr_ObjVal);
         gap = (ub-lb)/ub;
         cost_scenarios[s] = model[s]->get(GRB_DoubleAttr_ObjVal)/saainstance->getWeightTestScenarios(s);
      }
      std::vector<std::string> write() const {
         std::vector<std::string> output; 
         return output; 
      }
};


#endif
