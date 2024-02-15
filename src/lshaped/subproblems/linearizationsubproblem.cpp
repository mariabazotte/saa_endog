#include "linearizationsubproblem.hpp"

LinearizationSubproblem::LinearizationSubproblem(const AllSceInstance * inst, const Input & input) : 
                            AbstractSubProblem(input), instance(inst) {
    defineParameters();
    create();
}

LinearizationSubproblem::~LinearizationSubproblem(){
    delete[] recursive_one;
    delete[] recursive_two;
    delete[] mccormick_one;

    delete[] arc_capacity;
    delete[] flow_clients;
    delete[] flow_facilities;
    delete[] u_;
    delete[] v_;
    delete[] l_;

    for(int i = 0; i < count_mccormick; ++i){
        delete[] mccormick_two[i];
        delete[] mccormick_three[i];
        delete[] mccormick_four[i];
        delete[] dual_mccormick_three[i];
        delete[] dual_mccormick_four[i];
    }
    delete[] mccormick_two;
    delete[] mccormick_three;
    delete[] mccormick_four;
    delete[] dual_mccormick_three;
    delete[] dual_mccormick_four;

    delete[] z;
    delete[] w;
    for(int facility=0;facility<instance->getNbFacilities();++facility){
        for(int protec=0;protec<instance->getNbFacilityProtLevels();++protec){
            delete[] v[facility][protec];
        }
        delete[] q[facility];
        delete[] v[facility];
    }
    delete[] q;
    delete[] v;   
}

void LinearizationSubproblem::defineParameters() {
    AbstractSubProblem::defineParameters();
    model->set(GRB_IntParam_NumericFocus,1);
    model->set(GRB_DoubleParam_OptimalityTol,0.000000001);
}

void LinearizationSubproblem::create(){
    createVariables();
    createConstraints();
    initializeVectors();
}

void LinearizationSubproblem::initializeVectors(){
    u_ = new double[count_u];
    v_ = new double[count_v];
    l_ = new double[count_l];
    dual_mccormick_three = new double*[count_mccormick];
    dual_mccormick_four = new double*[count_mccormick];
    for(int i=0;i<count_mccormick;++i){
        dual_mccormick_three[i] = new double[count_mccormick_three];
        dual_mccormick_four[i] = new double[count_mccormick_four];
    }
}

void LinearizationSubproblem::createVariables(){
    int nb_facilities = instance->getNbFacilities();
    int nb_node_protec = instance->getNbFacilityProtLevels();

    z = model->addVars(instance->getNbArcsWithDummy(),GRB_CONTINUOUS);

    q = new GRBVar*[nb_facilities];
    v = new GRBVar**[nb_facilities];
    for(int f = 0; f < nb_facilities; ++f){
        q[f] = model->addVars(nb_node_protec,GRB_CONTINUOUS);
        v[f] = new GRBVar*[nb_node_protec];
        for(int p = 0; p < nb_node_protec; ++p){
            q[f][p].set(GRB_DoubleAttr_LB,0.0);
            q[f][p].set(GRB_StringAttr_VarName,"q[" + std::to_string(f) + "," + std::to_string(p) + "]");
            
            v[f][p] =  model->addVars(instance->getNbArcsWithDummy(),GRB_CONTINUOUS);
            for(int a = 0; a < instance->getNbArcsWithDummy(); ++a){
                v[f][p][a].set(GRB_DoubleAttr_LB,0.0);
                v[f][p][a].set(GRB_StringAttr_VarName,"v[" + std::to_string(f) + "," + std::to_string(p) + ","+ std::to_string(a) +"]");
            }
        }
    }

    for(int a=0;a<instance->getNbArcsWithDummy();++a){
        z[a].set(GRB_DoubleAttr_LB,0.0);
        z[a].set(GRB_StringAttr_VarName,"z[" + std::to_string(a) + "]");
    }

    for(int a=0;a<instance->getNbArcs();++a){
        z[a].set(GRB_DoubleAttr_UB,instance->getArcMaximumCapacity());
    }

    // Objective
    for(int p = 0; p < nb_node_protec; ++p){
        q[nb_facilities-1][p].set(GRB_DoubleAttr_Obj,1.0);
    }
}
        
void LinearizationSubproblem::createConstraints(){
    int nb_protec_node = instance->getNbFacilityProtLevels();
    int nb_facilities = instance->getNbFacilities();

    // Second-stage constraints
    recursive_one = new GRBConstr[nb_protec_node];
    recursive_two = new GRBConstr[nb_facilities-1];
    mccormick_one = new GRBConstr[nb_facilities*nb_protec_node];
    mccormick_two = new GRBConstr*[nb_facilities];
    
    count_mccormick = nb_facilities;
    count_mccormick_three = nb_protec_node*instance->getNbArcsWithDummy();
    count_mccormick_four = nb_protec_node*instance->getNbArcsWithDummy();

    mccormick_three = new GRBConstr*[count_mccormick];
    mccormick_four = new GRBConstr*[count_mccormick];
    for(int f = 0; f < nb_facilities; ++f){
        mccormick_two[f] = new GRBConstr[nb_protec_node*instance->getNbArcsWithDummy()];
        mccormick_three[f] = new GRBConstr[count_mccormick_three];
        mccormick_four[f] = new GRBConstr[count_mccormick_four];
    }

    count_u = instance->getNbClients();
    count_v = nb_facilities;
    count_l = instance->getNbArcs();

    arc_capacity = new GRBConstr[count_l];
    flow_clients = new GRBConstr[count_u];
    flow_facilities = new GRBConstr[count_v];

    // Clients demand
    for(int client=0;client<instance->getNbClients();++client){
        double demand = instance->getClientDemand(client);
        int n = instance->getNodeClient(client);
        GRBLinExpr flow = 0;
        for(int a: instance->getArcsArriving(n)){
            flow += z[a];
        }
        for(int a: instance->getArcsLeaving(n)){
            flow -= z[a];
        }
        flow_clients[client] = model->addConstr(flow == demand, "Flow_Client_" + std::to_string(client));
    }

    // Facilities capacity (do not consider dummy facility)
    for(int facility=0;facility<nb_facilities;++facility){
        int n = instance->getNodeFacility(facility);
        GRBLinExpr flow;
        for(int a: instance->getArcsArriving(n)){
            flow -= z[a];
        }
        for(int a: instance->getArcsLeaving(n)){
            flow += z[a];
        }
        flow_facilities[facility] = model->addConstr(flow <= 0, "Flow_Facility_" + std::to_string(facility));
    }

    // Arcs capacity (do not consider dummy facilities)
    for(int a=0;a<instance->getNbArcs();++a){
        GRBLinExpr capacity = z[a];
        arc_capacity[a] = model->addConstr(capacity <=0, "Capacity_Arc_" + std::to_string(a));
    }

    // Recursive constraints
    for(int p = 0; p < nb_protec_node; ++p){
        GRBLinExpr exp = q[0][p];
        for(int a = 0; a < instance->getNbArcsWithDummy(); ++a){
            exp += v[0][p][a];
        }
        recursive_one[p] = model->addConstr(exp == 0, "init_recursive_protec_" + std::to_string(p));
    }

    for(int f = 1; f < nb_facilities; ++f){
        GRBLinExpr exp = 0;
        for(int p = 0; p < nb_protec_node; ++p){
            exp += (q[f-1][p] - q[f][p]);
        }
        recursive_two[f-1] = model->addConstr(exp == 0, "recursive_facility_" + std::to_string(f));
    }

    // Linearization constraints
    int nb_arcs_dummy = instance->getNbArcsWithDummy();

    for(int f = 0; f < nb_facilities; ++f){
        for(int p = 0; p < nb_protec_node; ++p){
            GRBLinExpr exp = q[f][p];
            for(int a = 0; a < nb_arcs_dummy; ++a){
                exp += (-instance->getArcLength(a)*v[f][p][a]);
                GRBLinExpr exp2 = v[f][p][a] - z[a];
                GRBLinExpr exp3 = v[f][p][a];
                GRBLinExpr exp4 = v[f][p][a] - z[a];
                mccormick_two[f][p*nb_arcs_dummy+a] = model->addConstr(exp2 <= 0, "mccormick_2_facility" + std::to_string(f) +  "_protec_" + std::to_string(p) +  "_arc_" + std::to_string(a));
                mccormick_three[f][p*nb_arcs_dummy+a] = model->addConstr(exp3 <= 0, "mccormick_3_facility" + std::to_string(f) +  "_protec_" + std::to_string(p) +  "_arc_" + std::to_string(a));
                mccormick_four[f][p*nb_arcs_dummy+a] = model->addConstr(exp4 >= 0, "mccormick_4_facility" + std::to_string(f) +  "_protec_" + std::to_string(p) +  "_arc_" + std::to_string(a));
            }
            mccormick_one[f*nb_protec_node+p] = model->addConstr(exp <= 0, "mccormick_1_facility_" + std::to_string(f) +  "_protec_" + std::to_string(p));
        }
    }
}

void LinearizationSubproblem::solve(LinearizationMainProblem* main, double time, int scenario, bool ppd){
    int nb_facilities = instance->getNbFacilities();
    int nb_node_protec = instance->getNbFacilityProtLevels();
    for(int p = 0; p < nb_node_protec; ++p){
        q[nb_facilities-1][p].set(GRB_DoubleAttr_Obj,1.0);
    }
    
    for(int f = 0; f < instance->getNbFacilities(); ++f){
        flow_facilities[f].set(GRB_DoubleAttr_RHS,instance->getCapaFacility(scenario,f));
    }

    for(int p = 0; p < instance->getNbFacilityProtLevels(); ++p){
        for(int a = 0; a < instance->getNbArcsWithDummy(); ++a){
            model->chgCoeff(recursive_one[p], v[0][p][a], -(instance->getProbFacility(scenario,0,p)*instance->getArcLength(a)));
        }
    }

    for(int f = 1; f < instance->getNbFacilities(); ++f){
        for(int p = 0; p < instance->getNbFacilityProtLevels(); ++p){
            model->chgCoeff(recursive_two[f-1], q[f][p], -(1.0/std::max(0.00001,instance->getProbFacility(scenario,f,p))));
            //model->chgCoeff(recursive_two[f-1], q[f-1][p], instance->getProbFacility(scenario,f,p));
        }
    }

    if(ppd){
        for(int a=0;a<instance->getNbArcs();++a){
            int e = instance->getEdgeFromArc(a);
            arc_capacity[a].set(GRB_DoubleAttr_RHS,std::max(0.0,main->getY0(e)*instance->getArcMaximumCapacity()));
        }

        for(int facility = 0;facility<instance->getNbFacilities();++facility){
            for(int protec=0;protec<instance->getNbFacilityProtLevels();++protec){
                for(int a=0;a<instance->getNbArcsWithDummy();++a){
                    mccormick_three[facility][protec*instance->getNbArcsWithDummy()+a].set(GRB_DoubleAttr_RHS,std::max(0.0,main->getX0(facility,protec)*instance->getArcMaximumCapacity()));
                    mccormick_four[facility][protec*instance->getNbArcsWithDummy()+a].set(GRB_DoubleAttr_RHS,-(1.0-main->getX0(facility,protec))*instance->getArcMaximumCapacity());
                }
            }
        }
    }else{
        for(int a=0;a<instance->getNbArcs();++a){
            int e = instance->getEdgeFromArc(a);
            arc_capacity[a].set(GRB_DoubleAttr_RHS,std::max(0.0,main->getY(e)*instance->getArcMaximumCapacity()));
        }

        for(int f = 0; f < instance->getNbFacilities(); ++f){
            for(int p = 0; p < instance->getNbFacilityProtLevels(); ++p){
                for(int a = 0; a < instance->getNbArcsWithDummy(); ++a){
                    mccormick_three[f][p*instance->getNbArcsWithDummy()+a].set(GRB_DoubleAttr_RHS, std::max(0.0,main->getX(f,p)*instance->getArcMaximumCapacity()));
                    mccormick_four[f][p*instance->getNbArcsWithDummy()+a].set(GRB_DoubleAttr_RHS,-((1.0-main->getX(f,p))*instance->getArcMaximumCapacity()));
                }
            }
        }
    }
    model->set(GRB_DoubleParam_TimeLimit, time);
    model->update();
    model->optimize();

    if(input.getVerbose() == 3){
        model->write("linearizationsub" + std::to_string(scenario) + ".lp");
        model->write("linearizationsub" + std::to_string(scenario) + ".sol");
    }

    delete[] u_;
    delete[] v_;
    delete[] l_;

    for(int f = 0; f < instance->getNbFacilities(); ++f){
        delete[] dual_mccormick_three[f];
        delete[] dual_mccormick_four[f];
    }

    status =  statusFromGurobi(model->get(GRB_IntAttr_Status));
    if (status == Status::Unbounded || status == Status::Infeasible_or_Unbounded || status == Status::Infeasible){
        std::cerr << "ATTENTION: The Linearization sub problem should be always feasible. The status is " << status << ".\n";
        std::string file = input.getSolutionFile();
        if(status == Status::Infeasible){
            model->computeIIS();
            model->write(file.substr(0, file.size()-4) + "_linearizationsub.ilp"); 
        }
        model->write(file.substr(0, file.size()-4) + "_linearizationsub" + std::to_string(scenario) + ".lp");
        
        objective = GRB_INFINITY;
        u_ = model->get(GRB_DoubleAttr_UnbdRay,flow_clients,count_u);
        v_ = model->get(GRB_DoubleAttr_UnbdRay,flow_facilities,count_v);
        l_ = model->get(GRB_DoubleAttr_UnbdRay,arc_capacity,count_l);

        for(int i = 0; i < count_mccormick; ++i){
            dual_mccormick_three[i] = model->get(GRB_DoubleAttr_UnbdRay,mccormick_three[i],count_mccormick_three);
            dual_mccormick_four[i] = model->get(GRB_DoubleAttr_UnbdRay,mccormick_four[i],count_mccormick_four);
        }
        
        // May be infeasible because of problems with precision ???
        // throw std::runtime_error("ATTENTION: The Linearization sub problem should be always feasible.\n");
        // exit(0);    
    }
    else{
        u_ = model->get(GRB_DoubleAttr_Pi,flow_clients,count_u);
        v_ = model->get(GRB_DoubleAttr_Pi,flow_facilities,count_v);
        l_ = model->get(GRB_DoubleAttr_Pi,arc_capacity,count_l);

        for(int i = 0; i < count_mccormick; ++i){
            dual_mccormick_three[i] = model->get(GRB_DoubleAttr_Pi,mccormick_three[i],count_mccormick_three);
            dual_mccormick_four[i] = model->get(GRB_DoubleAttr_Pi,mccormick_four[i],count_mccormick_four);
        }
        objective = model->get(GRB_DoubleAttr_ObjVal)*instance->getProbDisruption(scenario);
    }
}
