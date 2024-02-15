#include "linearizationproblem.hpp"

LinearizationProblem::LinearizationProblem(const Input & input): 
                                       SimpleAbstractSolver(input,"Linearization"), 
                                       allinstance(AllSceInstance(input)){
   ev_lb = -GRB_INFINITY;
   ev_ub = GRB_INFINITY;
   ev_gap = 100.0;
   ev_time = 0.0;
   eev_lb = -GRB_INFINITY;
   eev_ub = GRB_INFINITY;
   eev_gap = 100.0;
   eev_time = 0.0;
   defineParameters();
   create();
}

LinearizationProblem::~LinearizationProblem(){
   int nb_facilities = allinstance.getNbFacilities();
   int nb_node_protec = allinstance.getNbFacilityProtLevels();
   int nb_arcs = allinstance.getNbArcsWithDummy();

   for(int f = 0; f < nb_facilities; ++f){
      for(int p = 0; p < nb_node_protec; ++p){
         for(int a = 0; a < nb_arcs; ++a){
            delete[] v[f][p][a];
         }
         delete[] q[f][p];
         delete[] v[f][p];
      }
      delete[] q[f];
      delete[] v[f];
   }
   delete[] q;
   delete[] v;
   delete[] x;
   delete[] y;
   delete[] w;
   delete[] z;
}

void LinearizationProblem::defineParameters(){
   SimpleAbstractSolver::defineParameters();
   model->set(GRB_IntParam_NumericFocus,1);
   model->set(GRB_DoubleParam_MIPGap,0.000000001);
   model->set(GRB_DoubleParam_OptimalityTol,0.000000001);
}

void LinearizationProblem::create(){
   createVariables();
   createConstraints();
}

void LinearizationProblem::createVariables(){
   int nb_facilities = allinstance.getNbFacilities();
   int nb_node_protec = allinstance.getNbFacilityProtLevels();
   int nb_sce = allinstance.getTotalNbScenarios();

   int count_x = nb_facilities*nb_node_protec;
   int count_y = allinstance.getNbEdges();
   int count_z = allinstance.getNbArcsWithDummy()*nb_sce;
   
   x = model->addVars(count_x,GRB_BINARY);
   y = model->addVars(count_y,GRB_BINARY);
   z = model->addVars(count_z,GRB_CONTINUOUS);

   q = new GRBVar**[nb_facilities];
   v = new GRBVar***[nb_facilities];
   for(int f = 0; f < nb_facilities; ++f){
      q[f] = new GRBVar*[nb_node_protec];
      v[f] = new GRBVar**[nb_node_protec];
      for(int p = 0; p < nb_node_protec; ++p){
         q[f][p] = model->addVars(nb_sce,GRB_CONTINUOUS);
         v[f][p] = new GRBVar*[allinstance.getNbArcsWithDummy()];
         for(int a = 0; a < allinstance.getNbArcsWithDummy(); ++a){
            v[f][p][a] = model->addVars(nb_sce,GRB_CONTINUOUS);
         }
      }
   }

   if(allinstance.getConstOrObjCost() == true){
      for(int f = 0; f < nb_facilities; ++f){
         for(int p = 0; p < nb_node_protec; ++p){
            x[f*nb_node_protec+p].set(GRB_DoubleAttr_Obj,allinstance.getCostFacilityProtection(f,p));
         }
      }
      for (int e = 0; e < allinstance.getNbEdges(); ++e){
         y[e].set(GRB_DoubleAttr_Obj,allinstance.getCostEdge(e));
      }
   } 

   for(int a = 0; a < allinstance.getNbArcs(); ++a){
      for(int s = 0; s < nb_sce; ++s){
         z[a*nb_sce+s].set(GRB_DoubleAttr_UB,allinstance.getArcMaximumCapacity());
      }
   }

   for(int p = 0; p < nb_node_protec; ++p){
      for(int s = 0; s < nb_sce; ++s){
         q[nb_facilities-1][p][s].set(GRB_DoubleAttr_Obj,allinstance.getProbDisruption(s));
      }
   }

   // Naming variables
   for (int f = 0; f < allinstance.getNbFacilities(); ++f){
      for(int p = 0; p < nb_node_protec; ++p){
         x[f*nb_node_protec+p].set(GRB_StringAttr_VarName,"x[" + std::to_string(f) + "," + std::to_string(p) + "]");
      }
   }

   for (int e = 0; e < allinstance.getNbEdges(); ++e){
      y[e].set(GRB_StringAttr_VarName,"y[" + std::to_string(e) + "]");
   }

   for(int a = 0; a < allinstance.getNbArcsWithDummy(); ++a){
      std::pair<int,int> arc = allinstance.getArc(a);
      for(int s = 0; s < nb_sce; ++s){
         z[a*nb_sce+s].set(GRB_StringAttr_VarName,"z[" + std::to_string(s) + "," + std::to_string(arc.first) + "," + std::to_string(arc.second) + "]");
      }
   }
}

void LinearizationProblem::createConstraints(){
   int nb_protec_node = allinstance.getNbFacilityProtLevels();
   int nb_facilities = allinstance.getNbFacilities();
   int nb_sce = allinstance.getTotalNbScenarios();

   /* First-stage constraints. */ 

   // Facility protection investment
   for(int f = 0; f < nb_facilities; ++f){
      GRBLinExpr invest = 0;
      for (int p = 0; p < nb_protec_node; ++p){
         invest += x[f*nb_protec_node+p];
      }
      model->addConstr(invest == 1, "investment_facility_" + std::to_string(f));
   }

   // Total cost/budget
   if(allinstance.getConstOrObjCost() == false){
      GRBLinExpr cost = 0;
      for(int f = 0; f < nb_facilities; ++f){
         for(int p = 0; p < nb_protec_node; ++p){
            cost += x[f*nb_protec_node+p]*allinstance.getCostFacilityProtection(f,p);
         }
      }
      for (int e = 0; e < allinstance.getNbEdges(); ++e){
         cost += y[e]*allinstance.getCostEdge(e);
      }
      model->addConstr(cost <= allinstance.getBudget(), "budget");
   }

   // Second-stage constraints

   // Clients demand
   for(int c = 0; c < allinstance.getNbClients(); ++c){
      double demand = allinstance.getClientDemand(c);
      int n = allinstance.getNodeClient(c);
      for(int s = 0; s < nb_sce; ++s){
         GRBLinExpr flow;
         for(int a: allinstance.getArcsArriving(n)){
            flow += z[a*nb_sce+s];
         }
         for(int a: allinstance.getArcsLeaving(n)){
            flow -= z[a*nb_sce+s];
         }
         model->addConstr(flow == demand, "flow_client_" + std::to_string(c) + "_scenario_" + std::to_string(s));
      }
   }

   // Facilities capacity (do not consider dummy facility)
   for(int f = 0; f < nb_facilities; ++f){
      int n = allinstance.getNodeFacility(f);
      for(int s = 0; s < nb_sce; ++s){
         GRBLinExpr flow;
         double capacity = allinstance.getCapaFacility(s,f);
         for(int a: allinstance.getArcsArriving(n)){
            flow -= z[a*nb_sce+s];
         }
         for(int a: allinstance.getArcsLeaving(n)){
            flow += z[a*nb_sce+s];
         }
         model->addConstr(flow <= capacity, "flow_facility_" + std::to_string(f) + "_scenario_" + std::to_string(s));
      }     
   }

   // Arcs capacity (do not consider dummy facilities)
   for(int a = 0; a < allinstance.getNbArcs(); ++a){
      for(int s = 0; s < nb_sce; ++s){
         int e = allinstance.getEdgeFromArc(a);
         GRBLinExpr capacity = z[a*nb_sce+s] - y[e]*allinstance.getArcMaximumCapacity();
         model->addConstr(capacity <= 0, "capacity_arc_" + std::to_string(a) + "_scenario_" + std::to_string(s));
      }
   }

   // Recursive constraints
   for(int p = 0; p < nb_protec_node; ++p){
      for(int s = 0; s < nb_sce; ++s){
         GRBLinExpr exp = q[0][p][s];
         for(int a = 0; a < allinstance.getNbArcsWithDummy(); ++a){
            //std::cout << allinstance.getArcLength(a) << std::endl;
            //std::cout << allinstance.getProbFacility(s,0,p) << std::endl;
            exp += (-allinstance.getProbFacility(s,0,p)*allinstance.getArcLength(a)*v[0][p][a][s]);
         }
         model->addConstr(exp == 0, "init_recursive_protec_" + std::to_string(p) + "_scenario_" + std::to_string(s));
      }
   }

   for(int f = 1; f < nb_facilities; ++f){
      for(int s = 0; s < nb_sce; ++s){
         GRBLinExpr exp  = 0;
         for(int p = 0; p < nb_protec_node; ++p){
            exp += (q[f-1][p][s] - (1.0/std::max(0.00001,allinstance.getProbFacility(s,f,p)))*q[f][p][s] );
            //exp += (allinstance.getProbFacility(s,f,p)*q[f-1][p][s] - q[f][p][s]);
         }
         model->addConstr(exp == 0, "recursive_facility_" + std::to_string(f) + "_scenario_" + std::to_string(s));
      }
   }

   // Linearization constraints
   for(int f = 0; f < nb_facilities; ++f){
      for(int p = 0; p < nb_protec_node; ++p){
         for(int s = 0; s < nb_sce; ++s){
            GRBLinExpr exp = q[f][p][s];
            for(int a = 0; a < allinstance.getNbArcsWithDummy(); ++a){
               exp += (-allinstance.getArcLength(a)*v[f][p][a][s]); 
            }
            model->addConstr(exp <= 0, "mccormick_1_facility_" + std::to_string(f) +  "_protec_" + std::to_string(p) + "_scenario_" + std::to_string(s));
         }
      }
   }

   for(int f = 0; f < nb_facilities; ++f){
      for(int p = 0; p < nb_protec_node; ++p){
         for(int a = 0; a < allinstance.getNbArcsWithDummy(); ++a){
            for(int s = 0; s < nb_sce; ++s){
               GRBLinExpr exp = v[f][p][a][s] - z[a*nb_sce+s];
               GRBLinExpr exp2 = v[f][p][a][s] - x[f*nb_protec_node+p]*allinstance.getArcMaximumCapacity();
               GRBLinExpr exp3 = v[f][p][a][s] - z[a*nb_sce+s] + (1.0 - x[f*nb_protec_node+p])*allinstance.getArcMaximumCapacity();

               model->addConstr(exp <= 0, "mccormick_2_facility_" + std::to_string(f) +  "_protec_" + std::to_string(p) +  "_arc_" + std::to_string(a) + "_scenario_" + std::to_string(s));
               model->addConstr(exp2 <= 0, "mccormick_3_facility_" + std::to_string(f) +  "_protec_" + std::to_string(p) +  "_arc_" + std::to_string(a) + "_scenario_" + std::to_string(s));
               model->addConstr(exp3 >= 0, "mccormick_4_facility_" + std::to_string(f) +  "_protec_" + std::to_string(p) +  "_arc_" + std::to_string(a) + "_scenario_" + std::to_string(s));
            }
         }
      }
   }
}

void LinearizationProblem::solve(){
   SimpleAbstractSolver::solve();
   solveEEV();
}

void LinearizationProblem::solveEEV(){
   int nb_facilities = allinstance.getNbFacilities();
   int nb_protec_node = allinstance.getNbFacilityProtLevels();

   GRBEnv *ev_env = new GRBEnv();
   GRBModel *ev_model = new GRBModel(*ev_env);
   
   GRBVar* ev_x = ev_model->addVars(nb_facilities*nb_protec_node,GRB_BINARY);
   GRBVar* ev_y = ev_model->addVars(allinstance.getNbEdges(),GRB_BINARY);
   GRBVar* ev_z =  ev_model->addVars(allinstance.getNbArcsWithDummy(),GRB_CONTINUOUS);

   // Objective
   for(int a = 0; a < allinstance.getNbArcsWithDummy(); ++a){
      ev_z[a].set(GRB_DoubleAttr_Obj,allinstance.getArcLength(a));
   }

   // Facility protection investment
   for(int f = 0; f < nb_facilities; ++f){
      GRBLinExpr invest = 0;
      for (int p = 0; p < nb_protec_node; ++p){
         invest += ev_x[f*nb_protec_node+p];
      }
      ev_model->addConstr(invest == 1, "investment_facility_" + std::to_string(f));
   }

   // Total cost/budget
   if(allinstance.getConstOrObjCost() == true){
      for(int f = 0; f < nb_facilities; ++f){
         for(int p = 0; p < nb_protec_node; ++p){
            ev_x[f*nb_protec_node+p].set(GRB_DoubleAttr_Obj,allinstance.getCostFacilityProtection(f,p));
         }
      }
      for (int e = 0; e < allinstance.getNbEdges(); ++e){
         ev_y[e].set(GRB_DoubleAttr_Obj,allinstance.getCostEdge(e));
      }
   } else {
      GRBLinExpr cost = 0;
      for(int f = 0; f < nb_facilities; ++f){
         for(int p = 0; p < nb_protec_node; ++p){
            cost += ev_x[f*nb_protec_node+p]*allinstance.getCostFacilityProtection(f,p);
         }
      }
      for (int e = 0; e < allinstance.getNbEdges(); ++e){
         cost += ev_y[e]*allinstance.getCostEdge(e);
      }
      ev_model->addConstr(cost <= allinstance.getBudget(), "budget");
   }

   for(int a = 0; a < allinstance.getNbArcs(); ++a){
      int e = allinstance.getEdgeFromArc(a);
      ev_model->addConstr(ev_z[a] <= ev_y[e]*allinstance.getArcMaximumCapacity(), "capacity_arc_" + std::to_string(a));
   }

   for(int c = 0; c < allinstance.getNbClients(); ++c){
      double demand = allinstance.getClientDemand(c);
      int n = allinstance.getNodeClient(c);
      GRBLinExpr flow = 0;
      for(int a: allinstance.getArcsArriving(n)){
         flow += ev_z[a];
      }
      for(int a: allinstance.getArcsLeaving(n)){
         flow -= ev_z[a];
      }
      ev_model->addConstr(flow == demand, "flow_client_" + std::to_string(c));
   }

   for(int f = 0; f < allinstance.getNbFacilities(); ++f){
      int n = allinstance.getNodeFacility(f);
      GRBLinExpr flow;
      for(int a: allinstance.getArcsArriving(n)){
         flow -= ev_z[a];
      }
      for(int a: allinstance.getArcsLeaving(n)){
         flow += ev_z[a];
      }
      for(int p = 0; p < nb_protec_node; ++p){
         flow -= ev_x[f*nb_protec_node+p]*allinstance.getMeanCapacity(f,p);
      }
      ev_model->addConstr(flow <= 0, "flow_facility_" + std::to_string(f));
   }
   ev_model->update();
   ev_model->optimize();

   ev_status = statusFromGurobi(ev_model->get(GRB_IntAttr_Status));
   if(ev_status == Status::Optimal){
      if(input.getVerbose() == 2){
         ev_model->write("evproblem.sol");
      }
      ev_lb = ev_model->get(GRB_DoubleAttr_ObjBound);
      ev_ub = ev_model->get(GRB_DoubleAttr_ObjVal);
      ev_gap = ev_model->get(GRB_DoubleAttr_MIPGap);
      ev_time = ev_model->get(GRB_DoubleAttr_Runtime);

      // Solving the EEV model
      double *ev_x_ = ev_model->get(GRB_DoubleAttr_X,ev_x,nb_facilities*nb_protec_node);
      double *ev_y_ = ev_model->get(GRB_DoubleAttr_X,ev_y,allinstance.getNbEdges());

      for(int i = 0; i < nb_facilities*nb_protec_node; ++i){
         model->addConstr(x[i] == ev_x_[i]);
      }
      for(int i = 0; i < allinstance.getNbEdges(); ++i){
         model->addConstr(y[i] == ev_y_[i]);
      }

      model->update();
      model->optimize();

      eev_status = statusFromGurobi(model->get(GRB_IntAttr_Status));
      if(eev_status == Status::Optimal){
         if(input.getVerbose() == 2){
            model->write("eevproblem.sol");
         }
         eev_lb = model->get(GRB_DoubleAttr_ObjBound);
         eev_ub = model->get(GRB_DoubleAttr_ObjVal);
         eev_gap = model->get(GRB_DoubleAttr_MIPGap);
         eev_time = model->get(GRB_DoubleAttr_Runtime);
      }
      delete[] ev_x_;
      delete[] ev_y_;
   }else{

   }

   delete[] ev_x;
   delete[] ev_y;
   delete[] ev_z;
   delete ev_model;
   delete ev_env; 
}