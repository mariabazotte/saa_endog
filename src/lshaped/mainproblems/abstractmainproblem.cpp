#include "abstractmainproblem.hpp"

AbstractMainProblem::AbstractMainProblem(const Input &input, int nb_scenarios) : 
                                        input(input), env(new GRBEnv()), model(new GRBModel(*env)),
                                        nb_subproblems(nb_scenarios), papadakos(input.getPapadakos()),
                                        callback_(input.getUseCallback()), user_cut(input.getUserAndLazy()),
                                        valid_inequality_EV(input.getUseValidInequalityEV()){
    lb = -GRB_INFINITY;
    ub = GRB_INFINITY;
    gap = 1.0;
    nb_branch_nodes = 0;
    status = Status::Not_Solved;
    defineParameters();
}

AbstractMainProblem::~AbstractMainProblem(){
    delete[] x_;
    delete[] y_;
    delete[] x0_;
    delete[] y0_;
    delete[] theta_;
    delete[] x;
    delete[] y;
    delete[] theta;
    delete model;
    delete env;
}

void AbstractMainProblem::initializeVectors(){
    x_ = new double[count_x];
    y_ = new double[count_y];
    x0_ = new double[count_x];
    y0_ = new double[count_y];
    theta_ = new double[count_theta];
}

void AbstractMainProblem::initializeCorePoint(){
    for (int i = 0; i < count_x; ++i){
        if(x_[i] <= 0.0 + intfeastol && x_[i] >= 0.0 - intfeastol)
            x0_[i] = 0;
        else if(x_[i] <= 1.0 + intfeastol && x_[i] >= 1.0 - intfeastol)
            x0_[i] = 1;
        else 
            x0_[i] = x_[i];
    }
    for (int i = 0; i < count_y; ++i){
        if(y_[i] <= 0.0 + intfeastol && y_[i] >= 0.0 - intfeastol)
            y0_[i] = 0;
        else if(y_[i] <= 1.0 + intfeastol && y_[i] >= 1.0 - intfeastol)
            y0_[i] = 1;
        else
            y0_[i] = y_[i];
    }
}

void AbstractMainProblem::updateCorePoint(){
    double lambda = 0.5;
    for (int i = 0; i < count_x; ++i){
        double x;
        if(x_[i] <= 0.0 + intfeastol && x_[i] >= 0.0 - intfeastol)
            x = 0;
        else if(x_[i] <= 1.0 + intfeastol && x_[i] >= 1.0 - intfeastol)
            x = 1;
        else
            x = x_[i];
        x0_[i] = x0_[i] * lambda + x * (1.00 - lambda);
    }
    for (int i = 0; i < count_y; ++i){
        double y;
        if(y_[i] <= 0.0 + intfeastol && y_[i] >= 0.0 - intfeastol)
            y = 0;
        else if(y_[i] <= 1.0 + intfeastol && y_[i] >= 1.0 - intfeastol)
            y = 1;
        else
            y = y_[i];
        y0_[i] = y0_[i] *lambda + y * (1.00 - lambda);
    }
}

void AbstractMainProblem::verifyVariables(){
    // Avoid problems with Integer Feasibility tolerance
    for(int i = 0; i < count_x; ++i){
        if(x_[i] <= 0.0 + intfeastol && x_[i] >= 0.0 - intfeastol)
            x_[i] = 0;
        else if(x_[i] <= 1.0 + intfeastol && x_[i] >= 1.0 - intfeastol)
            x_[i] = 1;
    }
    for(int i = 0; i < count_y; ++i){
        if(y_[i] <= 0.0 + intfeastol && y_[i] >= 0.0 - intfeastol)
            y_[i] = 0;
        else if(y_[i] <= 1.0 + intfeastol && y_[i] >= 1.0 - intfeastol)
            y_[i] = 1;
    }
}

void AbstractMainProblem::deleteVariables(){
    delete[] x_;
    delete[] y_; 
    delete[] theta_;
}

void AbstractMainProblem::defineParameters(){
    model->set(GRB_IntParam_Threads, input.getNbThreads());
    model->set(GRB_DoubleParam_TimeLimit, input.getTimeLimit());
    model->set(GRB_IntParam_IntegralityFocus,1); // Avoid problems with intfeastol
    model->set(GRB_IntParam_NumericFocus,1);
    model->set(GRB_DoubleParam_IntFeasTol,0.000000001);
    model->set(GRB_DoubleParam_MIPGap,0.0000001);
    intfeastol =  model->get(GRB_DoubleParam_IntFeasTol);
    if (callback_ == false) model->set(GRB_DoubleParam_MIPGap,0);
    if (input.getVerbose() <= 1) model->set(GRB_IntParam_OutputFlag, input.getVerbose());
    if (callback_ == true) model->set(GRB_IntParam_LazyConstraints, 1);
    if (callback_ == false && input.getVerbose() == 1) model->set(GRB_IntParam_OutputFlag, 0); 
}

void AbstractMainProblem::setLinearRelaxation(){
    for (int i = 0; i < count_x; ++i){
        x[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
        x[i].set(GRB_DoubleAttr_LB, 0.0);
        x[i].set(GRB_DoubleAttr_UB, 1.0);
    }
    for (int i = 0; i < count_y; ++i){
        y[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
        y[i].set(GRB_DoubleAttr_LB, 0.0);
        y[i].set(GRB_DoubleAttr_UB, 1.0);
    }
}

void AbstractMainProblem::setInteger(){
    for (int i = 0; i < count_x; ++i){
        x[i].set(GRB_CharAttr_VType, GRB_BINARY);
    }
    for (int i = 0; i < count_y; ++i){
        y[i].set(GRB_CharAttr_VType, GRB_BINARY);
    }
}

int AbstractMainProblem::getNbFoundSolutions(){
    return model->get(GRB_IntAttr_SolCount);
}

void AbstractMainProblem::solveCallback(){
    updateTimeLimit();

    model->setCallback(this);
    model->update();
    model->optimize();

    std::string output;
    status = statusFromGurobi(model->get(GRB_IntAttr_Status));
    if (status == Status::Optimal){
        lb = model->get(GRB_DoubleAttr_ObjBound);
        ub = model->get(GRB_DoubleAttr_ObjVal);
        gap = model->get(GRB_DoubleAttr_MIPGap);
        nb_branch_nodes = model->get(GRB_DoubleAttr_NodeCount);
        output = std::string("OPTIMAL: The optimal solution LShaped with callback is: ") + std::to_string(ub) +  std::string(".\n");
    }
    else if (status == Status::Time_Limit){
        lb = model->get(GRB_DoubleAttr_ObjBound);
        ub = model->get(GRB_DoubleAttr_ObjVal);
        gap = model->get(GRB_DoubleAttr_MIPGap);
        nb_branch_nodes = model->get(GRB_DoubleAttr_NodeCount);
        output = std::string("TIME: The best feasible solution found for the LShaped with callback is ") + std::to_string(ub) + std::string(".\n");
    }
    if (status == Status::Unbounded){
        output = std::string("UNBOUNDED: The problem LShaped with callback is unbounded.\n");
        std::cerr << "UNBOUNDED: The problem LShaped with callback is unbounded.\n";
    }
    if ((status == Status::Infeasible_or_Unbounded) && (status == Status::Infeasible)){
        output = std::string("GRB_INF_OR_UNBD: The problem LShaped with callback is unbounded or infeasible.\n");
        std::cerr << "GRB_INF_OR_UNBD: The problem LShaped with callback is unbounded or infeasible.\n";
    }
    if (input.getVerbose() >= 1)
        std::cout << output;
    if (input.getVerbose() == 2){
        model->write("main_problem.lp");
        model->write("main_problem.sol");
    }
}

void AbstractMainProblem::solve(){
    updateTimeLimit();
    model->update();
    model->optimize();

    std::string output;
    status = statusFromGurobi(model->get(GRB_IntAttr_Status));

    if (status == Status::Optimal){
        output = std::string("LShaped iteration: The optimal solution is: ") + std::to_string(model->get(GRB_DoubleAttr_ObjVal)) + std::string(".\n");
    }
    if (status == Status::Unbounded){
        output = std::string("LShaped iteration: The problem is unbounded.\n");
        throw std::runtime_error("LShaped iteration: The problem is unbounded.\n");
    }
    if ((status == Status::Infeasible) && (status == Status::Infeasible)){
        output = std::string("LShaped iteration: The problem is unbounded or infeasible.\n");
        throw std::runtime_error("LShaped iteration: The problem is unbounded or infeasible.\n");
    }
    if (input.getVerbose() == 2){
        std::cout << output;
        model->write("main_problem.lp");
        model->write("main_problem.sol");
    }
    if (status == Status::Time_Limit){
        lb = std::max(lb, model->get(GRB_DoubleAttr_ObjBound));
        return;
    }
    else{
        lb = model->get(GRB_DoubleAttr_ObjBound);
    }
    updateVariables();
    status = Status::Optimal;
}

void AbstractMainProblem::updateVariables(){
    deleteVariables();
    x_ = model->get(GRB_DoubleAttr_X,x, count_x);
    y_ = model->get(GRB_DoubleAttr_X,y, count_y);
    theta_ = model->get(GRB_DoubleAttr_X, theta, count_theta);
    sup = model->get(GRB_DoubleAttr_ObjVal) - std::accumulate(theta_, theta_+count_theta, 0.00000);
    verifyVariables();
}

void AbstractMainProblem::solveSubProblems(bool ppd){
    for (int s = 0; s < nb_subproblems; ++s){
        if(getTimeLimit() == 0.0){
            break;
        }
        sub->solve(this, getTimeLimit(), s, ppd);
        sup += sub->getObjective();
        if (sub->getStatus() != Status::Time_Limit && sub->getObjective() > theta_[s]){
            GRBLinExpr cut;
            separate(cut, s, ppd);
            if (sub->getObjective() >= GRB_INFINITY) {
                if (ppd)
                    model->addConstr(cut <= 0, "Viability_cut_ppd_" + std::to_string(nb_ppd_via_cuts));
                else
                    model->addConstr(cut <= 0, "Viability_cut_" + std::to_string(nb_via_cuts));
            }
            else {
                if (ppd)
                    model->addConstr(cut <= 0, "Optimality_cut_ppd_" + std::to_string(nb_ppd_opt_cuts));
                else
                    model->addConstr(cut <= 0, "Optimality_cut_" + std::to_string(nb_opt_cuts));
            }
        }
    }
}

void AbstractMainProblem::createLShapedVariables(){
    // L-shaped variables
    count_theta = nb_subproblems;
    theta = model->addVars(count_theta,GRB_CONTINUOUS);
    for(int i=0;i<count_theta;++i){
        theta[i].set(GRB_DoubleAttr_LB, 0.0);
        theta[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
        theta[i].set(GRB_DoubleAttr_Obj, 1);
        theta[i].set(GRB_StringAttr_VarName, "theta[" + std::to_string(i) + "]");
    }
    for(int i = 0; i < count_theta; ++i) ub_callback += theta[i];
}

void AbstractMainProblem::create(){
    createVariables();
    createConstraints();
    model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model->update();
    initializeVectors();
}

void AbstractMainProblem::createTrustRegion(){
    GRBLinExpr exp_x;
    GRBLinExpr exp_y;
    for(int i = 0; i < count_x; ++i){
        exp_x += x[i];
    }
    for(int i = 0; i < count_y; ++i){
        exp_y += y[i];
    }
    model->addConstr(exp_x <= 0, "trust_region_x");
    model->addConstr(exp_y <= 0, "trust_region_y");
    model->update();
}

void AbstractMainProblem::removeTrustRegion(){
    model->remove(model->getConstrByName("trust_region_x"));
    model->remove(model->getConstrByName("trust_region_y"));
    model->update();
}

void AbstractMainProblem::updateTrustRegion(){
    double delta = 0.4;
    int count_x_1 = 0;
    int count_y_1 = 0;
    for(int i = 0; i < count_x; ++i){
        if(x_[i] == 0)
            model->chgCoeff(model->getConstrByName("trust_region_x"), x[i], 1);
        if(x_[i] == 1){
            model->chgCoeff(model->getConstrByName("trust_region_x"), x[i], -1);
            count_x_1 += 1;
        }
        model->getConstrByName("trust_region_x").set(GRB_DoubleAttr_RHS, delta * count_x - count_x_1);
    }

    for(int i = 0; i < count_y; ++i){
        if(y_[i] == 0)
            model->chgCoeff(model->getConstrByName("trust_region_y"), y[i], 1);
        if(x_[i] == 1){
            model->chgCoeff(model->getConstrByName("trust_region_y"), y[i], -1);
            count_y_1 += 1;
        }
        model->getConstrByName("trust_region_y").set(GRB_DoubleAttr_RHS, delta * count_y - count_y_1);
    }
    model->update();
}

void AbstractMainProblem::separate(GRBLinExpr & cut, int s, bool ppd){
    if (ppd){
        if (sub->getObjective() >= GRB_INFINITY) {
            ++nb_ppd_via_cuts;
        }
        else{
            ++nb_ppd_opt_cuts;
            cut += (-theta[s]);
        }
    }
    else{
        if (sub->getObjective() >= GRB_INFINITY) {
            ++nb_via_cuts;
        }
        else{
            ++nb_opt_cuts;
            cut += (-theta[s]);
        }
    }
}

void AbstractMainProblem::callback(){
    try{
        if (where == GRB_CB_MIPSOL){   // Integer solution
            deleteVariables();
            x_ = getSolution(x, count_x);
            y_ = getSolution(y, count_y);
            theta_ = getSolution(theta, count_theta);
            verifyVariables();

            sup = getDoubleInfo(GRB_CB_MIPSOL_OBJ) - std::accumulate(theta_, theta_+count_theta, 0.00000);
            for (int s = 0; s < nb_subproblems; ++s){
                sub->solve(this, getTimeLimit(), s, false);
                sup += sub->getObjective();
                if (sub->getStatus() != Status::Time_Limit && sub->getObjective() > theta_[s]){
                    GRBLinExpr cut ;
                    separate(cut, s, false);
                    addLazy(cut <= 0);
                }   
                setSolution(theta[s],sub->getObjective());
            }
            setSolution(x,x_,count_x);
            setSolution(y,y_,count_y);
            useSolution();
            
            if(papadakos == true && not_init_ppd == true){
                initializeCorePoint();
                not_init_ppd = false;
            }else if(papadakos == true) {
                updateCorePoint();
                for (int s = 0; s < nb_subproblems; ++s){
                    sub->solve(this, getTimeLimit(), s, true);
                    if (sub->getStatus() != Status::Time_Limit && sub->getObjective() > theta_[s]){
                        GRBLinExpr cut ;
                        separate(cut, s, true);
                        addLazy(cut <= 0);
                    } 
                }
            }
            
            // addLazy(ub_callback <= sup); // Trying to use upper bound to prune some nodes
        }
        else if (where == GRB_CB_MIPNODE) {    // Fractional solution
            if (getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL && user_cut == true){
                if(getDoubleInfo(GRB_CB_MIPNODE_NODCNT) > 0) {
                    deleteVariables();
                    x_ = getNodeRel(x, count_x);
                    y_ = getNodeRel(y, count_y);
                    theta_ = getNodeRel(theta, count_theta);
                    for (int s = 0; s < nb_subproblems; ++s){
                        sub->solve(this, getTimeLimit(), s, false);
                        if(sub->getStatus() != Status::Time_Limit){
                            GRBLinExpr cut;
                            separate(cut, s, false);
                            addLazy(cut <= 0);
                        }
                    }
                }
            }
        }
    }
    catch (GRBException e){
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    catch (...){
        std::cout << "Error during callback" << std::endl;
    }
}

void AbstractMainProblem::updateTimeLimit(){
    model->set(GRB_DoubleParam_TimeLimit, std::max(0.0, input.getTimeLimit() - (time(NULL) - start_time)));
}

