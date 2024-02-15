#include "saatestproblem.hpp"

SAATest::SAATest(const Input & input,const SAAInstance* saainstance): AbstractSolver(input,"SAA Test"),
                                                                      saainstance(saainstance),instance(saainstance->getInstance()),
                                                                      num_scenarios(saainstance->getNbTestScenarios()),env(new GRBEnv()) { 
   model = new GRBModel*[num_scenarios];
   for(int i=0;i<num_scenarios;++i){
      model[i] = new GRBModel(*env);
   }
   z = new GRBVar*[num_scenarios];
   w = new GRBVar*[num_scenarios];
   cost_scenarios.resize(num_scenarios, 0.0);
   defineParameters();
   create(); 

   exact_problem = (input.getDistribution() != Input::Distribution::Normal && 
                   (instance->getNbDisrType()*std::pow(instance->getNbCapaStates(),instance->getNbFacilities())) < 20000);
   if(exact_problem){
      allsceinstance = new AllSceInstance(input);
      exact_env = new GRBEnv();
      exact_model = new GRBModel(*env);
      createExactTestProblem();  
   }
}

void SAATest::defineParameters(){
   for(int i=0;i<num_scenarios;++i){
      model[i]->set(GRB_IntParam_Threads, input.getNbThreads());
      model[i]->set(GRB_DoubleParam_TimeLimit, input.getTimeLimit());
      model[i]->set(GRB_DoubleParam_MIPGap, 1e-7);
      if(input.getVerbose()<=1) model[i]->set(GRB_IntParam_OutputFlag, input.getVerbose());
      model[i]->set(GRB_IntParam_OutputFlag, 0);
   }
}

void SAATest::solve(){
   lb = 0.0;
   ub = 0.0;
   for(int s = 0; s < num_scenarios; ++s){
      model[s]->update();
      model[s]->optimize();
      int gurobi_status = model[s]->get(GRB_IntAttr_Status);
      status = std::max(status,statusFromGurobi(gurobi_status));
      if((status == Status::Optimal) || (status == Status::Time_Limit)){
         updateBounds(s);
      }else{
         std::cout << "The SAA test problem is neither optimal nor time limit. Gurobi:" << gurobi_status << std::endl;
         if (status == Status::Unbounded){
            std::cerr << "The SAA problem is unbounded, this cannot happen according to the problem definition." << std::endl;
            throw std::runtime_error("SAA test problem is unbounded.\n");
            exit(0);
         }
         if ((status == Status::Infeasible_or_Unbounded) || (status == Status::Infeasible)){
            std::cerr << "The SAA test problem is unbounded or infeasible, this cannot happen according to the problem definition." << std::endl;
            model[s]->computeIIS();
            model[s]->write("saa_test" + std::to_string(s) + ".ilp");
            throw std::runtime_error("SAA problem is infeasible.\n");
            exit(0); 
         } 
      }
      if(input.getVerbose() >= 1){
         std::string output;
         if(status == Status::Optimal){
            output = "OPTIMAL: The optimal solution of the " + name_problem + " was " + std::to_string(ub) + ".\n";
         }else if(status == Status::Time_Limit){
            output = "TIME: The best feasible solution found for the " + name_problem + " was " + std::to_string(ub) + ".\n";
         }
         std::cout << output; 
         if(input.getVerbose() == 2){
            model[s]->write(name_problem + ".lp");
            model[s]->write(name_problem + ".sol");
         }
      }
   }
}

void SAATest::create(){
   for(int s = 0; s < num_scenarios; ++s){
      createProblem(s);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void SAATest::createProblem(int s){
   int count_z = instance->getNbArcsWithDummy();
   int count_w = instance->getNbFacilities();

   z[s] = model[s]->addVars(count_z,GRB_CONTINUOUS);
   w[s] = model[s]->addVars(count_w,GRB_CONTINUOUS);

   for(int a = 0; a < instance->getNbArcsWithDummy(); ++a){
      double length = instance->getArcLength(a);
      std::pair<int,int> arc = instance->getArc(a);
      z[s][a].set(GRB_DoubleAttr_Obj,length*saainstance->getWeightTestScenarios(s));
      z[s][a].set(GRB_StringAttr_VarName,"test_z[" + std::to_string(s) + "," + std::to_string(arc.first) + "," + std::to_string(arc.second) + "]");
   }

   for(int n = 0; n < instance->getNbFacilities(); ++n){
      w[s][n].set(GRB_DoubleAttr_LB,0.0);
   }

   for(int a = 0; a < instance->getNbArcs(); ++a){
      z[s][a].set(GRB_DoubleAttr_UB,instance->getArcMaximumCapacity());
   }

   // Clients demand 
   for(int client=0;client<instance->getNbClients();++client){
      double demand = instance->getClientDemand(client);
      int n = instance->getNodeClient(client);
      GRBLinExpr flow;
      for(int a: instance->getArcsArriving(n)){
         flow += z[s][a];
      }
      for(int a: instance->getArcsLeaving(n)){
         flow -= z[s][a];
      }
      model[s]->addConstr(flow == demand, "Flow_Client_" + std::to_string(client) + "_Scenario_" + std::to_string(s));
   }

   // Facilities capacity (do not consider dummy facility)
   for(int f = 0; f < instance->getNbFacilities(); ++f){
      int n = instance->getNodeFacility(f);
      GRBLinExpr flow;
      for(int a: instance->getArcsArriving(n)){
         flow -= z[s][a];
      }
      for(int a: instance->getArcsLeaving(n)){
         flow += z[s][a];
      }
      flow -= w[s][f];
      model[s]->addConstr(flow <= 0, "Flow_Facility_" + std::to_string(f) + 
                                       "_Scenario_" + std::to_string(s));     
   }
}

void SAATest::createExactTestProblem(){
   int nb_facilities = allsceinstance->getNbFacilities();
   int nb_scenarios = allsceinstance->getTotalNbScenarios();

   exact_model->set(GRB_DoubleParam_OptimalityTol,0.000000001);
   exact_model->set(GRB_IntParam_Threads, input.getNbThreads());
   if(input.getVerbose()<=1) exact_model->set(GRB_IntParam_OutputFlag, input.getVerbose());

   int count_exact_z = allsceinstance->getNbArcsWithDummy()*nb_scenarios;
   exact_z = exact_model->addVars(count_exact_z,GRB_CONTINUOUS);

   // Clients demand
   for(int c = 0; c < allsceinstance->getNbClients(); ++c){
      double demand = allsceinstance->getClientDemand(c);
      int n = allsceinstance->getNodeClient(c);
      for(int s = 0; s < nb_scenarios; ++s){
         GRBLinExpr flow;
         for(int a: allsceinstance->getArcsArriving(n)){
            flow += exact_z[a*nb_scenarios+s];
         }
         for(int a: allsceinstance->getArcsLeaving(n)){
            flow -= exact_z[a*nb_scenarios+s];
         }
         exact_model->addConstr(flow == demand, "flow_client_" + std::to_string(c) + "_scenario_" + std::to_string(s));
      }
   }

   // Facilities capacity (do not consider dummy facility)
   for(int f = 0; f < nb_facilities; ++f){
      int n = allsceinstance->getNodeFacility(f);
      for(int s = 0; s < nb_scenarios; ++s){
         GRBLinExpr flow;
         double capacity = allsceinstance->getCapaFacility(s,f);
         for(int a: allsceinstance->getArcsArriving(n)){
            flow -= exact_z[a*nb_scenarios+s];
         }
         for(int a: allsceinstance->getArcsLeaving(n)){
            flow += exact_z[a*nb_scenarios+s];
         }
         exact_model->addConstr(flow <= capacity, "flow_facility_" + std::to_string(f) + "_scenario_" + std::to_string(s));
      }     
   }

   exact_model->update();
}

void SAATest::solve(double * x_, double * y_, double * prob_, double * mean_, double * stddev_, double ** utility_){
   if(exact_problem){
      allsceinstance->computeScenarioProbability(x_,prob_,const_cast<const double **> (utility_));
      solveExactTestProblem(x_,y_);
   }else{
      int pr = saainstance->getTestProblem();
      double ** capacity_ = saainstance->getSolutionCapacityVector(pr,x_,prob_,mean_,stddev_,const_cast<const double **> (utility_));
      solve(x_,y_,capacity_);
      for(int s = 0; s < saainstance->getNbTestScenarios(); ++s) delete[] capacity_[s];
      delete[] capacity_;
   }
}

void SAATest::solve(double * x_,double * y_,double ** capacity_){
   int nb_scenarios = saainstance->getNbTestScenarios();
   int nb_protec_node = instance->getNbFacilityProtLevels();
   int nb_facilities = instance->getNbFacilities();

   for(int s = 0; s < num_scenarios; ++s){
      model[s]->set(GRB_DoubleAttr_UB,w[s],capacity_[s],nb_facilities);
   }
   
   for(int a = 0; a < instance->getNbArcs(); ++a){
      int e = instance->getEdgeFromArc(a);
      for(int s = 0; s < nb_scenarios; ++s){
         z[s][a].set(GRB_DoubleAttr_UB,y_[e]*instance->getArcMaximumCapacity());
      }
   }

   SAATest::solve();
   if(instance->getConstOrObjCost() == true){
      if(status == Status::Optimal || status == Status::Time_Limit){
         for(int f = 0; f < nb_facilities; ++f){
            for(int p = 0; p < nb_protec_node; ++p){
               lb += instance->getCostFacilityProtection(f,p)*x_[f*nb_protec_node+p]; 
               ub += instance->getCostFacilityProtection(f,p)*x_[f*nb_protec_node+p];
               for(int s = 0; s < nb_scenarios; ++s){
                  cost_scenarios[s] += instance->getCostFacilityProtection(f,p)*x_[f*nb_protec_node+p];
               }
            }
         }
         for (int e = 0; e < instance->getNbEdges(); ++e){
            lb += instance->getCostEdge(e)*y_[e];
            ub += instance->getCostEdge(e)*y_[e];
            for(int s = 0; s < nb_scenarios; ++s){
               cost_scenarios[s] += instance->getCostEdge(e)*y_[e];
            }
         }
      }
   }
}

void SAATest::solveExactTestProblem(double * x_,double * y_){
   // Compute probability scenarios and update objective
   int nb_protec_node = allsceinstance->getNbFacilityProtLevels();
   int nb_scenarios = allsceinstance->getTotalNbScenarios();
   for(int a = 0; a < allsceinstance->getNbArcsWithDummy(); ++a){
      for(int s = 0; s < nb_scenarios; ++s){
         exact_z[a*nb_scenarios+s].set(GRB_DoubleAttr_Obj,allsceinstance->getProbabilityScenario(s)*allsceinstance->getArcLength(a));
      }
   }

   // Arcs capacity (do not consider dummy facilities)
   for(int a = 0; a < allsceinstance->getNbArcs(); ++a){
      for(int s = 0; s < nb_scenarios; ++s){
         int e = allsceinstance->getEdgeFromArc(a);
         exact_z[a*nb_scenarios+s].set(GRB_DoubleAttr_UB,y_[e]*allsceinstance->getArcMaximumCapacity());
      }
   }

   exact_model->update();
   if(input.getVerbose()>=1) exact_model->write("exactproblem.lp");
   exact_model->optimize();

   status = statusFromGurobi(exact_model->get(GRB_IntAttr_Status));
   if(status == Status::Optimal){
      if(input.getVerbose() == 2){
         exact_model->write("exactproblem.sol");
      }
      
      ub = exact_model->get(GRB_DoubleAttr_ObjVal);
      lb = exact_model->get(GRB_DoubleAttr_ObjBound);
      if(allsceinstance->getConstOrObjCost() == true){
         for (int f = 0; f < allsceinstance->getNbFacilities(); ++f){
            for(int p = 0; p < nb_protec_node; ++p){
               ub += x_[f*nb_protec_node+p]*allsceinstance->getCostFacilityProtection(f,p);
               lb += x_[f*nb_protec_node+p]*allsceinstance->getCostFacilityProtection(f,p);
            }
         }
         for (int e = 0; e < allsceinstance->getNbEdges(); ++e){
            ub += y_[e]*allsceinstance->getCostEdge(e);
            lb += y_[e]*allsceinstance->getCostEdge(e);
         }
      }
   }else{
      std:: cout << "Status" << status << std::endl;
      std::cerr << "The exact test problem is not optimal, this should not happen." << std::endl;
      throw std::runtime_error("Exact test problem is not optimal.\n");
      exit(0);
   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

