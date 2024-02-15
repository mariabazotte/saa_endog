#ifndef ABSTRACT_SOLVER_HPP
#define ABSTRACT_SOLVER_HPP

#include "gurobi_c++.h"
#include "../input/input.hpp"
#include "status.hpp"
#include "../data/saainstance.hpp"
#include "../data/instance.hpp"
#include "../data/allsceinstance.hpp"

class AbstractSolver {
    protected:
        const Input & input;
        std::string name_problem;

        double lb;
        double ub;
        double gap;
        double nb_branch_nodes;
        double time_;
        int int_status;
        Status status;

        virtual void create() = 0;
        virtual void defineParameters() = 0;

    public:
        AbstractSolver(const Input & input, std::string name) : input(input), name_problem(name) {
            lb = -GRB_INFINITY;
            ub = GRB_INFINITY;
            gap = 1.0;
            nb_branch_nodes = 0;
            time_ = 0.0;
            status = Status::Not_Solved;
        }

        virtual ~AbstractSolver(){}

        double getLB() const { return lb; }
        double getUB() const { return ub; }
        double getGap() const { return gap; }
        double getTime() const { return time_; }
        double getNbBranchNodes() const { return nb_branch_nodes; }
        Status getStatus() const { return status; }
        void setTime(double t) { time_ = t; }  

        virtual void solve() = 0;

        virtual std::vector<std::string> write() const = 0;

        std::string writeline() const {
            return std::to_string(status) + std::string(";") + 
                   std::to_string(lb) + std::string(";") + 
                   std::to_string(ub) + std::string(";") + 
                   std::to_string(gap) + std::string(";") + 
                   std::to_string(time_) + std::string(";") +
                   std::to_string(nb_branch_nodes) + std::string(";");  
        }
};

class SimpleAbstractSolver: public AbstractSolver {
    protected:
        GRBEnv *env = NULL;
        GRBModel *model = NULL;

        virtual void defineParameters(){
            model->set(GRB_IntParam_Threads, input.getNbThreads());
            model->set(GRB_DoubleParam_TimeLimit, input.getTimeLimit());
            if(input.getVerbose()<=1) model->set(GRB_IntParam_OutputFlag, input.getVerbose());
            if(input.getVerbose()<=1) std::cout << input.getVerbose() << std::endl;
        }

    public:
        SimpleAbstractSolver(const Input & input, std::string name) : AbstractSolver(input,name),
                                                    env(new GRBEnv()), model(new GRBModel(*env)) {}

        virtual ~SimpleAbstractSolver(){ delete model;
                                         delete env; }
        
        virtual void updateBounds(){
            lb = model->get(GRB_DoubleAttr_ObjBound);
            ub = model->get(GRB_DoubleAttr_ObjVal);
            gap = model->get(GRB_DoubleAttr_MIPGap);
            nb_branch_nodes = model->get(GRB_DoubleAttr_NodeCount);
        }

        virtual void solve() {
            model->update();
            model->optimize();
            int gurobi_status = model->get(GRB_IntAttr_Status);
            status = statusFromGurobi(gurobi_status);
            time_ = model->get(GRB_DoubleAttr_Runtime);
            std::string output;
            if(status == Status::Optimal){
                updateBounds();
                output = "OPTIMAL: The optimal solution of the " + name_problem + " was " + std::to_string(ub) + ".\n";
            }else if(status == Status::Time_Limit){
                updateBounds();
                output = "TIME: The best feasible solution found for the " + name_problem + " was " + std::to_string(ub) + ".\n";
            }
            if (status == Status::Unbounded){
                output = "UNBOUNDED: The " + name_problem + " is unbounded.\n";
            }
            if (status == Status::Infeasible){
                output = "INFEASIBLE: The " + name_problem + " is infeasible.\n";
                model->computeIIS();
                model->write(name_problem + ".ilp");
            }
            if ((status == Status::Infeasible_or_Unbounded) || (status == Status::Infeasible)){
                output = "INF_OR_UNBD: The " + name_problem + " is unbounded or infeasible.\n";
                model->computeIIS();
                model->write(name_problem + ".ilp");
            }
            if(input.getVerbose() >= 1)
                std::cout << output; 
            if(input.getVerbose() == 2){
                model->write(name_problem + ".lp");
                model->write(name_problem + ".sol");
            }
        }
};

#endif