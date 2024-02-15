#include "expectedproblem.hpp"

void ExpectedValueProblem::solve(){
    SimpleAbstractSolver::solve();
    getVariables();
    solveSAAEEV();
    if(input.getDistribution() != Input::Distribution::Normal && nb_endogenous_scenarios < 100000) solveEEV();
}

void ExpectedValueProblem::getVariables(){
    int nb_disr_type = instance->getNbDisrType();
    int capacity_levels = instance->getNbCapaStates();
    int nb_protec_node = instance->getNbFacilityProtLevels();
    int nb_facilities = instance->getNbFacilities();

    x_ = model->get(GRB_DoubleAttr_X,x,count_x);
    y_ = model->get(GRB_DoubleAttr_X,y,count_y);
    integerFeasibility();

    if(prob) prob_ = model->get(GRB_DoubleAttr_X,prob,nb_facilities*nb_disr_type);
    if(mean) mean_ = model->get(GRB_DoubleAttr_X,mean,nb_facilities*nb_disr_type);
    if(stddev) stddev_ = model->get(GRB_DoubleAttr_X,stddev,nb_facilities*nb_disr_type);
    if(aux_expected){
        utility_ = new double*[nb_facilities*nb_disr_type];
        sum_utility_ = new double[nb_facilities*nb_disr_type];
        for(int dt = 0; dt < nb_disr_type; ++dt){
            for(int n = 0; n < nb_facilities; ++n){
                sum_utility_[dt*nb_facilities+n] = 0.0;
                utility_[dt*nb_facilities+n] = new double[capacity_levels];
                for(int w = 0; w < capacity_levels; ++w){
                    utility_[dt*nb_facilities+n][w] = 0.0;
                    for(int nn = 0; nn < nb_facilities; ++nn){
                        double weight = (nn == n) ? 1.0 : instance->getFacilitiesImpact();
                        for(int p = 0; p < nb_protec_node; ++p){
                            utility_[dt*nb_facilities+n][w] += weight*instance->getUtility(dt,nn,p,w)*x_[nn*nb_protec_node+p];
                        }
                    }
                    sum_utility_[dt*nb_facilities+n] += utility_[dt*nb_facilities+n][w];
                }
                for(int w = 0; w < capacity_levels; ++w){
                    utility_[dt*nb_facilities+n][w] = utility_[dt*nb_facilities+n][w]/sum_utility_[dt*nb_facilities+n];
                }
            }
        }   
    }
}

void ExpectedValueProblem::solveSAAEEV(){
    int pr = saainstance->getTestProblem();
    double ** capacity_ = saainstance->getSolutionCapacityVector(pr,x_,prob_,mean_,stddev_,const_cast<const double **> (utility_));

    // Solving SAA test problem
    saatest->solve(x_,y_,capacity_);
    saa_eev_lb = saatest->getLB();
    saa_eev_ub = saatest->getUB();
    saa_eev_gap = saatest->getGap();
    saa_eev_time = saatest->getTime();
    saa_eev_status = saatest->getStatus();
    saa_eev_var_ub = 0.0;
    for(int s = 0; s < saainstance->getNbTestScenarios(); ++s){
        saa_eev_var_ub += std::pow((saatest->getCostScenario(s)-saa_eev_ub),2)*saainstance->getNbRepeatedTestScenarios(s);
    }
    long long N = saainstance->getOriginalNbTestScenarios();
    long long aux = (N*(N-1));
    if(N>1)saa_eev_var_ub = saa_eev_var_ub/aux;
    saa_eev_stat_ub = saa_eev_ub + saainstance->getCriticalNormal()*std::sqrt(saa_eev_var_ub);

    for(int s = 0; s < saainstance->getNbTestScenarios(); ++s) delete[] capacity_[s];
    delete[] capacity_;
}

void ExpectedValueProblem::integerFeasibility(){
    // Avoiding problems with integer feasibility tolerance
    double intfeastol = model->get(GRB_DoubleParam_IntFeasTol);
    for(int j = 0; j < count_x; ++j){
        if(x_[j] <= 0.0 + intfeastol && x_[j] >= 0.0 - intfeastol) x_[j] = 0;
        else if(x_[j] <= 1.0 + intfeastol && x_[j] >= 1.0 - intfeastol) x_[j] = 1;
    }        
    for(int j = 0; j < count_y; ++j){
        if(y_[j] <= 0.0 + intfeastol && y_[j] >= 0.0 - intfeastol) y_[j] = 0;
        else if(y_[j] <= 1.0 + intfeastol && y_[j] >= 1.0 - intfeastol) y_[j] = 1;
    }
}

void ExpectedValueProblem::createDistribution(){
    if(input.getDistribution() == Input::Distribution::DiscreteChoice)
        discreteChoiceDistribution();
    else if(input.getDistribution() == Input::Distribution::Normal)
        normalDistribution();
    else if(input.getDistribution() == Input::Distribution::Binomial)
        binomialDistribution();
    else if(input.getDistribution() == Input::Distribution::StdNormalization)
        stdNormaDistribution();
}

void ExpectedValueProblem::discreteChoiceDistribution(){
    int nb_facilities = instance->getNbFacilities();
    int nb_protec_node = instance->getNbFacilityProtLevels();

    // Facilities capacity (do not consider dummy facility)
    for(int f = 0; f < nb_facilities; ++f){
        int n = instance->getNodeFacility(f); 
        GRBLinExpr flow;
        for(int a: instance->getArcsArriving(n)){
            flow -= z[a];
        }
        for(int a: instance->getArcsLeaving(n)){
            flow += z[a];
        }
        for(int p = 0; p < nb_protec_node; ++p){
            flow -= x[f*nb_protec_node+p]*instance->getMeanCapacity(f,p);
        }
        model->addConstr(flow <= 0, "Flow_Facility_" + std::to_string(f));
    }
}

void ExpectedValueProblem::normalDistribution(){
    int nb_disr_type = instance->getNbDisrType();
    int nb_protec_node = instance->getNbFacilityProtLevels();
    int nb_facilities = instance->getNbFacilities();
    
    mean = model->addVars(nb_facilities*nb_disr_type,GRB_CONTINUOUS);
    stddev = model->addVars(nb_facilities*nb_disr_type,GRB_CONTINUOUS);

    // Limite values for mean and naming variables
    for(int dt = 0; dt < nb_disr_type; ++dt) {
        for(int n = 0; n < nb_facilities; ++n){
            mean[dt*nb_facilities+n].set(GRB_DoubleAttr_LB,0.0);

            mean[dt*nb_facilities+n].set(GRB_StringAttr_VarName,
                                    "mean[" + std::to_string(dt) + 
                                    "," + std::to_string(n) + "]");
            stddev[dt*nb_facilities+n].set(GRB_DoubleAttr_LB,0.0);

            stddev[dt*nb_facilities+n].set(GRB_StringAttr_VarName,
                                    "stddev[" + std::to_string(dt) + 
                                    "," + std::to_string(n) + "]");
        }
    }

    // Mean depend on the protection decisions
    double rho = 1 + (instance->getNbFacilities()-1)*instance->getFacilitiesImpact();
    for(int dt = 0; dt < nb_disr_type; ++dt) {
        for(int n = 0; n < nb_facilities; ++n){
            GRBLinExpr def_mean = rho*mean[dt*nb_facilities+n];
            GRBLinExpr def_stddev = rho*stddev[dt*nb_facilities+n];
            for(int nn = 0; nn < nb_facilities; ++nn){
                double weight = (nn == n) ? 1.0 : instance->getFacilitiesImpact();
                def_stddev -= weight*instance->getStdDevNormalDist(dt,nn);
                for(int p = 0; p < nb_protec_node; ++p){
                    def_mean -= (weight*instance->getMeanNormalDist(dt,nn,p)*x[nn*nb_protec_node+p]);
                }
            }
            model->addConstr(def_mean == 0.0 , "Mean_Def_Disr_" + std::to_string(dt) + "_Facility/Node_" + std::to_string(n));
            model->addConstr(def_stddev == 0.0 , "StdDev_Def_Disr_" + std::to_string(dt) + "_Facility/Node_" + std::to_string(n));
        }
    }

    // Facilities capacity (do not consider dummy facility)
    for(int f = 0; f < nb_facilities; ++f){
        int n = instance->getNodeFacility(f); 
        GRBLinExpr flow;
        for(int a: instance->getArcsArriving(n)){
            flow -= z[a];
        }
        for(int a: instance->getArcsLeaving(n)){
            flow += z[a];
        }
        for(int dt = 0; dt < nb_disr_type; ++dt){
            flow -= mean[dt*nb_facilities+f]*instance->getProbDisrType(dt);
        }
        if(instance->getUncertainNoDisr() == false){
            flow -= instance->getFullCapacity()*instance->getProbNoDisr();
        }
        model->addConstr(flow <= 0, "Flow_Facility_" + std::to_string(f));
    }
}

void ExpectedValueProblem::binomialDistribution(){
    int nb_disr_type = instance->getNbDisrType();
    int nb_facilities = instance->getNbFacilities();
    int nb_protec_node = instance->getNbFacilityProtLevels();
    
    prob = model->addVars(nb_facilities*nb_disr_type,GRB_CONTINUOUS);
    for(int dt = 0; dt < nb_disr_type; ++dt){
        for(int n = 0; n < nb_facilities; ++n){
            prob[dt*nb_facilities+n].set(GRB_DoubleAttr_LB,0.0);
            prob[dt*nb_facilities+n].set(GRB_DoubleAttr_UB,1.0);
        }
    }

    // Probability equal to weighted average
    double rho = 1 + (nb_facilities-1)*instance->getFacilitiesImpact();
    for(int dt = 0; dt < nb_disr_type; ++dt){
        for(int n = 0; n < nb_facilities; ++n){
            GRBLinExpr exp_prob = rho*prob[dt*nb_facilities+n];
            for(int nn = 0; nn < nb_facilities; ++nn){
                double weight = (nn == n) ? 1.0 : instance->getFacilitiesImpact();
                for(int p = 0; p < nb_protec_node; ++p){
                    exp_prob -= weight*instance->getBernoulliProbCapa(dt,nn,p)*x[nn*nb_protec_node+p];
                }
            }
            model->addConstr(exp_prob == 0, "Prob_Def_Disr_" + std::to_string(dt) + "_Facility_" + std::to_string(n));
        }
    }

    // Facilities capacity (do not consider dummy facility)
    for(int f = 0; f < nb_facilities; ++f){
        int n = instance->getNodeFacility(f); 
        GRBLinExpr flow;
        for(int a: instance->getArcsArriving(n)){
            flow -= z[a];
        }
        for(int a: instance->getArcsLeaving(n)){
            flow += z[a];
        }
        for(int dt = 0; dt < nb_disr_type; ++dt){
            flow -= prob[dt*nb_facilities+f]*instance->getFullCapacity()*instance->getProbDisrType(dt);
        }
        if(instance->getUncertainNoDisr() == false){
            flow -= instance->getFullCapacity()*instance->getProbNoDisr();
        }
        model->addConstr(flow <= 0, "Flow_Facility_" + std::to_string(f));
    }
}

void ExpectedValueProblem::stdNormaDistribution(){
    int nb_disr_type = instance->getNbDisrType();
    int nb_capa_levels = instance->getNbCapaStates();
    int nb_protec_node = instance->getNbFacilityProtLevels();
    int nb_facilities = instance->getNbFacilities();

    // Defining variables
    aux_expected = model->addVars(nb_facilities*nb_disr_type,GRB_CONTINUOUS);

    count_A = nb_facilities*nb_disr_type;

    A = new GRBVar*[nb_facilities*nb_disr_type];
    for(int dt = 0; dt < nb_disr_type; ++dt){
        for(int n = 0; n < nb_facilities; ++n){
            A[dt*nb_facilities+n] = model->addVars(nb_facilities*nb_protec_node,GRB_CONTINUOUS);
        }
    }

    // McCormick constraints
    for(int dt = 0; dt < nb_disr_type; ++dt){
        for(int n = 0; n < nb_facilities; ++n){
            for(int n2 = 0; n2 < nb_facilities; ++n2){
                for(int p = 0; p < nb_protec_node; ++p){
                    model->addConstr(A[dt*nb_facilities+n][n2*nb_protec_node+p] <= aux_expected[dt*nb_facilities+n], "McCormick1_Facility/Node_" + std::to_string(n) + 
                                                    "_Facility/Node_" + std::to_string(n2) + "_Protec_" + std::to_string(p));
                    model->addConstr(A[dt*nb_facilities+n][n2*nb_protec_node+p] <= instance->getFullCapacity()*x[n2*nb_protec_node+p], "McCormick2_Facility/Node_" + std::to_string(n) + 
                                                    "_Facility/Node_" + std::to_string(n2) + "_Protec_" + std::to_string(p));
                    model->addConstr(A[dt*nb_facilities+n][n2*nb_protec_node+p] >=  aux_expected[dt*nb_facilities+n] - instance->getFullCapacity()*(1-x[n2*nb_protec_node+p]), "McCormick2_Facility/Node_" + std::to_string(n) + 
                                                    "_Facility/Node_" + std::to_string(n2) + "_Protec_" + std::to_string(p));
                }
            }
        }
    }

    // Definition of the expected value for each node and disruption type
    for(int dt = 0; dt < nb_disr_type; ++dt){
        for(int n = 0; n < nb_facilities; ++n){
            GRBLinExpr expected_def;
            for(int w = 0; w < nb_capa_levels; ++w){
                int capacity = instance->getCapacity(w);
                for(int p = 0; p < nb_protec_node; ++p){
                    expected_def += instance->getUtility(dt,n,p,w)*A[dt*nb_facilities+n][n*nb_protec_node+p];
                    expected_def -= instance->getUtility(dt,n,p,w)*capacity*x[n*nb_protec_node+p];
                }
                for(int n2 = 0; n2 < nb_facilities; ++n2){
                    if(n2 != n){
                        for(int p = 0; p < nb_protec_node; ++p){
                            expected_def += instance->getUtility(dt,n2,p,w)*instance->getFacilitiesImpact()
                                            *A[dt*nb_facilities+n][n2*nb_protec_node+p];
                            expected_def -= instance->getUtility(dt,n2,p,w)*instance->getFacilitiesImpact()
                                            *capacity*x[n2*nb_protec_node+p];
                        }
                    }
                }
            }
            model->addConstr(expected_def == 0, "Expected_Val_Def_Disr_" + std::to_string(dt) + 
                                                "_Node/Facility_" + std::to_string(n));
        }
    }

    // Facilities capacity (do not consider dummy facility)
    for(int f = 0; f < nb_facilities; ++f){
        int n = instance->getNodeFacility(f); 
        GRBLinExpr flow;
        for(int a: instance->getArcsArriving(n)){
            flow -= z[a];
        }
        for(int a: instance->getArcsLeaving(n)){
            flow += z[a];
        }
        for(int dt = 0; dt < nb_disr_type; ++dt){
            flow -= aux_expected[dt*nb_facilities+f]*instance->getProbDisrType(dt);
        }
        if(instance->getUncertainNoDisr() == false){
            flow -= instance->getFullCapacity()*instance->getProbNoDisr();
        }
        model->addConstr(flow <= 0, "Flow_Facility_" + std::to_string(f));
    }
}

void ExpectedValueProblem::create(){
    int nb_facilities = instance->getNbFacilities();
    int nb_protec_node = instance->getNbFacilityProtLevels();

    count_x = nb_facilities*nb_protec_node;
    count_y = instance->getNbEdges();
    int count_z = instance->getNbArcsWithDummy();

    x = model->addVars(count_x,GRB_BINARY);
    y = model->addVars(count_y,GRB_BINARY);
    z = model->addVars(count_z,GRB_CONTINUOUS);

    // Objective
    for(int a = 0; a < instance->getNbArcsWithDummy(); ++a){
      double length = instance->getArcLength(a);
      std::pair<int,int> arc = instance->getArc(a);
        z[a].set(GRB_DoubleAttr_Obj, length);
        z[a].set(GRB_StringAttr_VarName,"z[" + std::to_string(arc.first) + 
                                        "," + std::to_string(arc.second) + "]");
    }

    if(instance->getConstOrObjCost() == true){
        for (int f = 0; f < instance->getNbFacilities(); ++f){
            for(int p = 0; p < nb_protec_node; ++p){
                x[f*nb_protec_node+p].set(GRB_DoubleAttr_Obj,instance->getCostFacilityProtection(f,p));
            }
        }
        for (int e = 0; e < instance->getNbEdges(); ++e){
            y[e].set(GRB_DoubleAttr_Obj,instance->getCostEdge(e));
        }
    } 

    // Facilities protection
    for (int f = 0; f < nb_facilities; ++f){
        GRBLinExpr invest = 0;
        for (int p = 0; p < nb_protec_node; ++p){
            invest += x[f*nb_protec_node+p];
        }
        model->addConstr(invest == 1, "Investment_Facility_" + std::to_string(f));
    }

    // Budget
    if(instance->getConstOrObjCost() == false){
        GRBLinExpr cost = 0;
        for (int f = 0; f < nb_facilities; ++f){
            for(int p = 0; p < nb_protec_node; ++p){
                cost += x[f*nb_protec_node+p]*instance->getCostFacilityProtection(f,p);
            }
        }
        for (int e = 0; e < instance->getNbEdges(); ++e){
            cost += y[e]*instance->getCostEdge(e);
        }
        model->addConstr(cost <= instance->getBudget(), "Budget");
    }

    // Clients demand 
    for(int c = 0; c < instance->getNbClients(); ++c){
        double demand = instance->getClientDemand(c);
        int n = instance->getNodeClient(c);
        GRBLinExpr flow;
        for(int a: instance->getArcsArriving(n)){
            flow += z[a];
        }
        for(int a: instance->getArcsLeaving(n)){
            flow -= z[a];
        }
        model->addConstr(flow == demand, "Flow_Client_" + std::to_string(c));
    }

    // Arcs capacity (do not consider dummy facilities)
    for(int a = 0; a < instance->getNbArcs(); ++a){
        int e = instance->getEdgeFromArc(a);
        GRBLinExpr capacity = z[a] - instance->getArcMaximumCapacity()*y[e];
        model->addConstr(capacity <= 0, "Capacity_Arc_" + std::to_string(a));
    }
    createDistribution();
}

void ExpectedValueProblem::solveEEV(){
    AllSceInstance allinstance = AllSceInstance(input);
    int nb_protec_node = allinstance.getNbFacilityProtLevels();
    int nb_facilities = allinstance.getNbFacilities();
    int nb_sce = allinstance.getTotalNbScenarios();
    eev_total_nb_scenarios = allinstance.getTotalNbScenarios();

    eev_model->set(GRB_DoubleParam_OptimalityTol,0.000000001);

    int count_eev_z = instance->getNbArcsWithDummy()*nb_sce;
    eev_z = eev_model->addVars(count_eev_z,GRB_CONTINUOUS);

    // Objective function
    // Compute probability scenarios 
    allinstance.computeScenarioProbability(x_,prob_,const_cast<const double **> (utility_));
    for(int a = 0; a < allinstance.getNbArcsWithDummy(); ++a){
        for(int s = 0; s < nb_sce; ++s){
            eev_z[a*nb_sce+s].set(GRB_DoubleAttr_Obj,allinstance.getProbabilityScenario(s)*allinstance.getArcLength(a));
        }
    }

    // Clients demand
    for(int c = 0; c < allinstance.getNbClients(); ++c){
        double demand = allinstance.getClientDemand(c);
        int n = allinstance.getNodeClient(c);
        for(int s = 0; s < nb_sce; ++s){
            GRBLinExpr flow;
            for(int a: allinstance.getArcsArriving(n)){
                flow += eev_z[a*nb_sce+s];
            }
            for(int a: allinstance.getArcsLeaving(n)){
                flow -= eev_z[a*nb_sce+s];
            }
            eev_model->addConstr(flow == demand, "flow_client_" + std::to_string(c) + "_scenario_" + std::to_string(s));
        }
    }

    // Facilities capacity (do not consider dummy facility)
    for(int f = 0; f < nb_facilities; ++f){
        int n = allinstance.getNodeFacility(f);
        for(int s = 0; s < nb_sce; ++s){
            GRBLinExpr flow;
            double capacity = allinstance.getCapaFacility(s,f);
            for(int a: allinstance.getArcsArriving(n)){
                flow -= eev_z[a*nb_sce+s];
            }
            for(int a: allinstance.getArcsLeaving(n)){
                flow += eev_z[a*nb_sce+s];
            }
            eev_model->addConstr(flow <= capacity, "flow_facility_" + std::to_string(f) + "_scenario_" + std::to_string(s));
        }     
    }

    // Arcs capacity (do not consider dummy facilities)
    for(int a = 0; a < allinstance.getNbArcs(); ++a){
        for(int s = 0; s < nb_sce; ++s){
            int e = allinstance.getEdgeFromArc(a);
            GRBLinExpr capacity = eev_z[a*nb_sce+s] - y_[e]*allinstance.getArcMaximumCapacity();
            eev_model->addConstr(capacity <= 0, "capacity_arc_" + std::to_string(a) + "_scenario_" + std::to_string(s));
        }
    }

    eev_model->update();
    eev_model->optimize();

    eev_status = statusFromGurobi(eev_model->get(GRB_IntAttr_Status));
    if(eev_status == Status::Optimal){
        if(input.getVerbose() == 2){
            eev_model->write("eevproblem.sol");
        }
        eev_lb = eev_model->get(GRB_DoubleAttr_ObjBound);
        eev_ub = eev_model->get(GRB_DoubleAttr_ObjVal);
        eev_gap = (eev_ub-eev_lb)/eev_ub;
        eev_time = eev_model->get(GRB_DoubleAttr_Runtime);

        if(allinstance.getConstOrObjCost() == true){
            for (int f = 0; f < allinstance.getNbFacilities(); ++f){
                for(int p = 0; p < nb_protec_node; ++p){
                    eev_lb += x_[f*nb_protec_node+p]*allinstance.getCostFacilityProtection(f,p);
                    eev_ub += x_[f*nb_protec_node+p]*allinstance.getCostFacilityProtection(f,p);
                }
            }
            for (int e = 0; e < allinstance.getNbEdges(); ++e){
                eev_lb += y_[e]*allinstance.getCostEdge(e);
                eev_ub += y_[e]*allinstance.getCostEdge(e);
            }
        }
    }
}
