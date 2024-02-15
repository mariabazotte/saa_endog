#ifndef SOLVER_FABRIC_H
#define SOLVER_FABRIC_H

#include "abstractsolver.hpp"
#include "saa/saaproblem.hpp"
#include "baseline/linearizationproblem.hpp"
#include "baseline/expectedproblem.hpp"
#include "../lshaped/lshaped.hpp"

class SolverFactory{
    public:
        inline AbstractSolver* createSolver(const Input & input){
            Input::Solver chosenSolver = input.getSolver();
            switch (chosenSolver)
            {
                case Input::Solver::SAA :{
                    return new SAAProblem(input);
                    break;
                }
                
                case Input::Solver::Linearization :{
                    return new LinearizationProblem(input);
                    break;         
                }
                
                case Input::Solver::Linearization_LShaped :{
                    return new LinearizationLShaped(input);
                    break;
                }

                case Input::Solver::Expected :{
                    return new ExpectedValueProblem(input);
                    break;
                }

                default:{
                    throw std::string("ERROR: Invalid Solver option in parameters file.");
                    exit(0);
                    break;
                }
            }
            return NULL;
        }
};

#endif