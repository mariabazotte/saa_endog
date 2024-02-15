#ifndef LINEAR_MAIN_PROBLEM_HPP
#define LINEAR_MAIN_PROBLEM_HPP

#include "abstractmainproblem.hpp"
#include "../../data/allsceinstance.hpp"
#include "../subproblems/linearizationsubproblem.hpp"

class LinearizationSubproblem;

class LinearizationMainProblem : public AbstractMainProblem{
    private:
        const AllSceInstance *allinstance;
        LinearizationSubproblem* lsub;

    protected:
        void separate(GRBLinExpr & cut,int s,bool ppd);
        void createVariables();
        void createConstraints();

    public:
        LinearizationMainProblem(const AllSceInstance *, const Input &);
        ~LinearizationMainProblem();

        double getY(int e) const { 
            if(y_[e] <= 0.0 + intfeastol && y_[e] >= 0.0 - intfeastol)
                return 0.0;
            else if(y_[e] <= 1.0 + intfeastol && y_[e] >= 1.0 - intfeastol)
                return 1.0;
            else{
                if(y_[e] < 0){
                    std::cout << y_[e] << std::endl;
                    return y_[e];
                }
                return y_[e]; 
            }
        }
        
        double getX(int n, int p) const { 
            if(x_[n*allinstance->getNbFacilityProtLevels()+p] <= 0.0 + intfeastol &&
               x_[n*allinstance->getNbFacilityProtLevels()+p] >= 0.0 - intfeastol ){
                return 0.0; 
            }else if(x_[n*allinstance->getNbFacilityProtLevels()+p] <= 1.0 + intfeastol &&
                    x_[n*allinstance->getNbFacilityProtLevels()+p] >= 1.0 - intfeastol ){
                return 1.0;
            }else{
                return x_[n*allinstance->getNbFacilityProtLevels()+p];
            }
        }
        
        double getY0(int e) const { return y0_[e]; }
        double getX0(int n, int p) const { return x0_[n*allinstance->getNbFacilityProtLevels()+p]; }     
};

#endif