#include "linearizationmainproblem.hpp"

LinearizationMainProblem::LinearizationMainProblem(const AllSceInstance *inst, const Input & input) : 
                                                AbstractMainProblem(input, inst->getTotalNbScenarios()),
                                                allinstance(inst){
    lsub = new LinearizationSubproblem(inst,input);                                             
    sub = lsub;
    create();  
}

LinearizationMainProblem::~LinearizationMainProblem() { 
    delete lsub; 
}  

void LinearizationMainProblem::createVariables(){
    AbstractMainProblem::createLShapedVariables();
    int nb_protec_node = allinstance->getNbFacilityProtLevels();
    int nb_facilities = allinstance->getNbFacilities();
    int nb_edges = allinstance->getNbEdges();

    count_x = nb_facilities*nb_protec_node;
    count_y = nb_edges;
    
    x = model->addVars(count_x,GRB_BINARY);
    y = model->addVars(count_y,GRB_BINARY);

    // Objective
    if(allinstance->getConstOrObjCost() == true){
        for (int f = 0; f < nb_facilities; ++f){
            for(int p = 0; p < nb_protec_node; ++p){
                x[f*nb_protec_node+p].set(GRB_DoubleAttr_Obj,allinstance->getCostFacilityProtection(f,p));
            }
        }
        for (int e = 0; e < nb_edges; ++e){
            y[e].set(GRB_DoubleAttr_Obj,allinstance->getCostEdge(e));
        }
    }   

    // Naming variables
    for (int f = 0; f < nb_facilities; ++f){
        for(int p = 0; p < nb_protec_node; ++p){
            x[f*nb_protec_node+p].set(GRB_StringAttr_VarName,"x[" + std::to_string(f) + "," + std::to_string(p) + "]");
        }
    }
    for (int e = 0; e < nb_edges; ++e){
        y[e].set(GRB_StringAttr_VarName,"y[" + std::to_string(e) + "]");
    }
}

void LinearizationMainProblem::createConstraints(){
    int nb_protec_node = allinstance->getNbFacilityProtLevels();
    int nb_facilities = allinstance->getNbFacilities();

    for (int f = 0; f < nb_facilities; ++f){
        GRBLinExpr invest;
        for (int p = 0; p < nb_protec_node; ++p){
            invest += x[f*nb_protec_node+p];
        }
        model->addConstr(invest == 1, "investment_facility_" + std::to_string(f) );
    }

    if(allinstance->getConstOrObjCost() == false){
        GRBLinExpr cost = 0;
        for (int f = 0; f < allinstance->getNbFacilities(); ++f){
            for(int p = 0; p < nb_protec_node; ++p){
                cost += x[f*nb_protec_node+p]*allinstance->getCostFacilityProtection(f,p);
            }
        }
        for (int edge = 0; edge < allinstance->getNbEdges(); ++edge){
            cost += y[edge]*allinstance->getCostEdge(edge);
        }
        model->addConstr(cost <= allinstance->getBudget(), "Budget");
    }

    if(input.getUseValidInequalityEV()){
        GRBVar* mean_z = model->addVars(allinstance->getNbArcsWithDummy(),GRB_CONTINUOUS);
        for(int a = 0; a < allinstance->getNbArcs(); ++a){
            int e = allinstance->getEdgeFromArc(a);
            model->addConstr(mean_z[a] <= y[e]*allinstance->getArcMaximumCapacity(), "mean_arc_capacity_" + std::to_string(a));
        }
        for(int c = 0; c < allinstance->getNbClients(); ++c){
            double demand = allinstance->getClientDemand(c);
            int n = allinstance->getNodeClient(c);
            GRBLinExpr flow = 0;
            for(int a: allinstance->getArcsArriving(n)){
                flow += mean_z[a];
            }
            for(int a: allinstance->getArcsLeaving(n)){
                flow -= mean_z[a];
            }
            model->addConstr(flow == demand, "mean_flow_client_" + std::to_string(c));
        }
        for(int f = 0; f < allinstance->getNbFacilities(); ++f){
            int n = allinstance->getNodeFacility(f);
            GRBLinExpr flow;
            for(int a: allinstance->getArcsArriving(n)){
                flow -= mean_z[a];
            }
            for(int a: allinstance->getArcsLeaving(n)){
                flow += mean_z[a];
            }
            for(int p = 0; p < nb_protec_node; ++p){
                flow -= x[f*nb_protec_node+p]*allinstance->getMeanCapacity(f,p);
            }
            model->addConstr(flow <= 0, "flow_facility_" + std::to_string(f));
        }
        GRBLinExpr lb = 0;
        for(int i = 0; i < count_theta; ++i){
            lb += theta[i];
        }
        for(int a = 0; a < allinstance->getNbArcsWithDummy(); ++a){
            lb -= allinstance->getArcLength(a)*mean_z[a];
        }
        model->addConstr(lb >= 0, "ev_cut");
        delete[] mean_z;
    }

    model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model->update();
}

void LinearizationMainProblem::separate(GRBLinExpr & cut,int s,bool ppd){
    int nb_node_protec = allinstance->getNbFacilityProtLevels();
    for(int f = 0; f < allinstance->getNbFacilities(); ++f){
        for(int p = 0; p < nb_node_protec; ++p){
            for(int a = 0; a < allinstance->getNbArcsWithDummy(); ++a){
                double mc_three = lsub->getDualMccormickThree(f,p,a);
                double mc_four = lsub->getDualMccormickFour(f,p,a);
                cut += allinstance->getProbDisruption(s)*x[f*nb_node_protec+p]*mc_three*allinstance->getArcMaximumCapacity();
                cut += allinstance->getProbDisruption(s)*x[f*nb_node_protec+p]*mc_four*allinstance->getArcMaximumCapacity();
                cut -= allinstance->getProbDisruption(s)*mc_four*allinstance->getArcMaximumCapacity();
            }
        }
    }
    for(int a = 0; a < allinstance->getNbArcs(); ++a){
        double l_ = lsub->getl_(a);
        int e = allinstance->getEdgeFromArc(a);
        cut += allinstance->getProbDisruption(s)*y[e]*l_*allinstance->getArcMaximumCapacity();
    }
    cut += allinstance->getProbDisruption(s)*allinstance->getRhsSimpleCutClient(lsub->getu_vector());
    cut += allinstance->getProbDisruption(s)*allinstance->getRhsSimpleCutFacility(s,lsub->getv_vector());
    AbstractMainProblem::separate(cut,s,ppd);
}

