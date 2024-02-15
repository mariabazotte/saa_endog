#include "saaproblem.hpp"

SAAProblem::SAAProblem(const Input & input) : AbstractSolver(input,"SAAProblem"), saainstance(new SAAInstance(input)), 
                                              instance(saainstance->getInstance()), saatest(new SAATest(input,saainstance)), env(new GRBEnv()) {
   model = new GRBModel*[saainstance->getNbProblems()];
   for(int i=0;i<saainstance->getNbProblems();++i){
      model[i] = new GRBModel(*env);
   }
   x = new GRBVar*[saainstance->getNbProblems()];
   y = new GRBVar*[saainstance->getNbProblems()];
   z = new GRBVar*[saainstance->getNbProblems()];

   if(input.getDistribution() == Input::Distribution::Binomial){
      prob = new GRBVar*[saainstance->getNbProblems()];
      capa = new GRBVar**[saainstance->getNbProblems()];
   }else if(input.getDistribution() == Input::Distribution::Normal){
      mean = new GRBVar*[saainstance->getNbProblems()];
      stddev = new GRBVar*[saainstance->getNbProblems()];
   }else if(input.getDistribution() == Input::Distribution::StdNormalization){
      capa = new GRBVar**[saainstance->getNbProblems()];
      cum_capa = new GRBVar**[saainstance->getNbProblems()];
      utility = new GRBVar**[saainstance->getNbProblems()];
      sum_utility = new GRBVar*[saainstance->getNbProblems()];
      mccormick1 = new GRBVar**[saainstance->getNbProblems()];
   }

   variance_lower = 0.0;
   variance_upper = 0.0;
   variance_gap = 0.0;
   statistical_lower = 0.0;
   statistical_upper = 0.0;
   statistical_gap = 0.0;
   var_nb_branch_nodes = 0.0;
   useinequality = input.getUseInequality();
   nb_saapr_not_sol_opt = 0;
   exact_problem = (input.getDistribution() != Input::Distribution::Normal && 
                   (instance->getNbDisrType()*std::pow(instance->getNbCapaStates(),instance->getNbFacilities())) < 20000);
   defineParameters();
   create();
}

SAAProblem::~SAAProblem(){ 
   for(int i = 0; i < saainstance->getNbProblems(); ++i) {   
      if(capa) {
         for(int s = 0; s < saainstance->getNbScenarios(i); ++s){
            delete[] capa[i][s];
         }
         delete[] capa[i];
      }
      if(cum_capa){
         for(int s = 0; s < saainstance->getNbScenarios(i); ++s){
            delete[] cum_capa[i][s];
         }
         delete[] cum_capa[i];
      }
      if(mccormick1){
         for(int s = 0; s < saainstance->getNbScenarios(i); ++s){
            delete[] mccormick1[i][s];
         }
         delete[] mccormick1[i];
      }
      if(utility){
         int m = instance->getNbDisrType()*instance->getNbFacilities();
         for(int j = 0; j < m; ++j){
            delete[] utility[i][j];
         }
         delete[] utility[i];
      }
      delete[] x[i]; 
      delete[] y[i]; 
      delete[] z[i];
      if(sum_utility) delete[] sum_utility[i];
      if(prob) delete[] prob[i];
      if(mean) delete[] mean[i];
      if(stddev) delete[] stddev[i];
      delete model[i]; 
   }
   delete[] x; 
   delete[] y; 
   delete[] z;
   if(mean) delete[] mean;
   if(stddev) delete[] stddev;
   if(prob) delete[] prob;
   if(capa) delete[] capa;
   if(cum_capa) delete[] cum_capa;
   if(utility) delete[] utility;
   if(mccormick1) delete[] mccormick1;
   if(sum_utility) delete[] sum_utility;
   delete[] model; 
   delete env;
   delete saatest; 
   delete saainstance;
}

void SAAProblem::computeEstimators(){
    variance_lower = 0.0;
    var_nb_branch_nodes = 0.0;

    long long M = saainstance->getNbProblems();
    lb = std::accumulate(lower_estimates.begin(),lower_estimates.end(),0.0)/M;
    nb_branch_nodes = std::accumulate(count_bnbnodes.begin(),count_bnbnodes.end(),0.0)/M;
    for(int i = 0; i < saainstance->getNbProblems(); ++i) {
        variance_lower += std::pow((lower_estimates[i]-lb),2);
        var_nb_branch_nodes += std::pow((count_bnbnodes[i]-nb_branch_nodes),2);
    }
    if(M > 1) variance_lower = variance_lower/(M*(M-1));
    if(M > 1) var_nb_branch_nodes = var_nb_branch_nodes/(M*(M-1));
    variance_gap = variance_lower + variance_upper;
    gap = (ub - lb)/ub;  

    // Statistical estimators considering 95% accuracy (alpha = 5 %)
    statistical_gap = (ub - lb)  + saainstance->getCriticalNormal()*std::sqrt(variance_gap);             
    statistical_lower = lb - saainstance->getCriticalTStudent()*std::sqrt(variance_lower);     
    statistical_upper = ub + saainstance->getCriticalNormal()*std::sqrt(variance_upper);
}

void SAAProblem::solve(){
    lb = 0.0;
    ub = GRB_INFINITY;
    lower_estimates.resize(saainstance->getNbProblems());
    count_bnbnodes.resize(saainstance->getNbProblems());
   start_time = time(NULL); 
    for(int i = 0; i < saainstance->getNbProblems(); ++i){
        Status current_status = solveProblem(i);
        updateIteration(current_status,i);
        printInfo(current_status,i);
        if(getNbFoundSolutions(i) > 0){
            updateVariables(i);
            solveTestProblem();
            deleteVariables();
        }
    }
    computeEstimators();
    time_ = time(NULL) - start_time;
}

void SAAProblem::solveTestProblem(){
   // Solving SAA test problem
   saatest->solve(x_,y_,prob_,mean_,stddev_,utility_);

   if(saatest->getStatus() == Status::Optimal && saatest->getUB() < ub){
      ub = std::min(ub, saatest->getUB());
      variance_upper = 0.0;
      if(exact_problem == false){
         for(int s = 0; s < saainstance->getNbTestScenarios(); ++s){
            variance_upper += std::pow((saatest->getCostScenario(s)-ub),2)*saainstance->getNbRepeatedTestScenarios(s);
         }
         long long N = saainstance->getOriginalNbTestScenarios();
         if(N > 1) variance_upper = variance_upper/(N*(N-1));
      }
   }
}

int SAAProblem::getNbFoundSolutions(int i){
   return model[i]->get(GRB_IntAttr_SolCount);
}

Status SAAProblem::solveProblem(int i){
   model[i]->set(GRB_DoubleParam_TimeLimit, std::max(0.0,input.getTimeLimit()-(time(NULL)-start_time)));
   model[i]->optimize();
   Status current_status = statusFromGurobi(model[i]->get(GRB_IntAttr_Status));
   return current_status;
}

void SAAProblem::updateIteration(Status current_status, int i){
   status = std::max(status,current_status);
   if(current_status == Status::Optimal || current_status == Status::Time_Limit){
      lower_estimates[i] = model[i]->get(GRB_DoubleAttr_ObjBound); 
      count_bnbnodes[i] = model[i]->get(GRB_DoubleAttr_NodeCount);
      if(current_status == Status::Time_Limit) ++nb_saapr_not_sol_opt;
   }else{
      if (current_status == Status::Unbounded){
         std::cerr << "The SAA problem is unbounded, this cannot happen according to the problem definition." << std::endl;
         throw std::runtime_error("SAA problem is unbounded.\n");
         exit(0); 
      }
      if ((current_status == Status::Infeasible_or_Unbounded) || (current_status == Status::Infeasible)){
         std::cerr << "The SAA problem is unbounded or infeasible, this cannot happen according to the problem definition." << std::endl;
         model[i]->computeIIS();
         model[i]->write("saa" + std::to_string(i) + ".ilp");
         throw std::runtime_error("SAA problem is infeasible.\n");
         exit(0); 
      }
   }
}

void SAAProblem::printInfo(Status current_status, int i){
   if(input.getVerbose() >= 1){
      std::string output;
      if(current_status == Status::Optimal){
         output = std::string("The optimal solution from problem ") + std::to_string(i) 
               + std::string(" is ") + std::to_string(model[i]->get(GRB_DoubleAttr_ObjVal)) + std::string(".\n");
      }else if(current_status == Status::Time_Limit){
         output = std::string("TIME: The best bound for problem ") +  std::to_string(i) 
               + std::string(" is ") + std::to_string(model[i]->get(GRB_DoubleAttr_ObjBound)) + std::string(".\n");
         output += std::string("TIME: The best feasible solution found for problem ") + std::to_string(i) 
               + std::string(" is ") + std::to_string(model[i]->get(GRB_DoubleAttr_ObjVal)) + std::string(".\n");
      }
      std::cout << output;
      if(input.getVerbose() == 2)
         model[i]->write("saa" + std::to_string(i) + ".sol");
   }
}

void SAAProblem::defineParameters(){
   for(int i=0;i<saainstance->getNbProblems();++i){
      model[i]->set(GRB_IntParam_Threads, input.getNbThreads());
      model[i]->set(GRB_DoubleParam_TimeLimit, input.getTimeLimit());
      model[i]->set(GRB_DoubleParam_MIPGap, 1e-7);
      // model[i]->set(GRB_DoubleParam_IntFeasTol, 1e-9);
      // model[i]->set(GRB_DoubleParam_FeasibilityTol, 1e-9);
      if(input.getVerbose()<=1) model[i]->set(GRB_IntParam_OutputFlag, input.getVerbose());
   }
   feasibility_tol =  model[0]->get(GRB_DoubleParam_FeasibilityTol);
}

void SAAProblem::create(){
   for(int i = 0; i < saainstance->getNbProblems(); ++i){
      generateCommonVariables(i);
      generateCommonConstraints(i);
      generateDistribution(i);
      model[i]->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
      model[i]->update();
      if(input.getVerbose() == 2){
         model[i]->write("saa_" + std::to_string(i)+ ".lp");
      }
   }  
}

void SAAProblem::generateDistribution(int i){
   if(saainstance->getDistribution() == Input::Distribution::DiscreteChoice) discreteChoiceDistribution(i);
   else if(saainstance->getDistribution() == Input::Distribution::Binomial) binomialDistribution(i);
   else if(saainstance->getDistribution() == Input::Distribution::Normal) normalDistribution(i);
   else if(saainstance->getDistribution() == Input::Distribution::StdNormalization) stdNormaDistribution(i);
}

void SAAProblem::updateVariables(int i){
   int nb_disr_type = instance->getNbDisrType();
   int capacity_levels = instance->getNbCapaStates();
   int nb_facilities = instance->getNbFacilities();

   x_ = model[i]->get(GRB_DoubleAttr_X,x[i], count_x);
   y_ = model[i]->get(GRB_DoubleAttr_X,y[i], count_y);

   if(prob) prob_ = model[i]->get(GRB_DoubleAttr_X,prob[i],nb_facilities*nb_disr_type);
   if(mean) mean_ = model[i]->get(GRB_DoubleAttr_X,mean[i],nb_facilities*nb_disr_type);
   if(stddev) stddev_ = model[i]->get(GRB_DoubleAttr_X,stddev[i],nb_facilities*nb_disr_type);
   if(sum_utility) sum_utility_ = model[i]->get(GRB_DoubleAttr_X,sum_utility[i],nb_facilities*nb_disr_type);
   if(utility) {
      utility_ = new double*[nb_facilities*nb_disr_type];
      for(int dt = 0; dt < nb_disr_type; ++dt){
         for(int n = 0; n < nb_facilities; ++n){
            utility_[dt*nb_facilities+n] = model[i]->get(GRB_DoubleAttr_X,utility[i][dt*nb_facilities+n],capacity_levels);
            for(int w = 0; w < capacity_levels; ++w){
               if(utility_[dt*nb_facilities+n][w] > -feasibility_tol && utility_[dt*nb_facilities+n][w] < 0.0)  utility_[dt*nb_facilities+n][w] = 0.0;
               utility_[dt*nb_facilities+n][w] = utility_[dt*nb_facilities+n][w]/sum_utility_[dt*nb_facilities+n];
            }
         } 
      }
   }

   // Avoiding problems with integer feasibility tolerance
   double intfeastol = model[i]->get(GRB_DoubleParam_IntFeasTol);
   for(int j = 0; j < count_x; ++j){
      if(x_[j] <= 0.0 + intfeastol && x_[j] >= 0.0 - intfeastol) x_[j] = 0;
      else if(x_[j] <= 1.0 + intfeastol && x_[j] >= 1.0 - intfeastol) x_[j] = 1;
   }        
   for(int j = 0; j < count_y; ++j){
      if(y_[j] <= 0.0 + intfeastol && y_[j] >= 0.0 - intfeastol) y_[j] = 0;
      else if(y_[j] <= 1.0 + intfeastol && y_[j] >= 1.0 - intfeastol) y_[j] = 1;
   }
}

void SAAProblem::deleteVariables(){
   delete[] x_;
   delete[] y_;
   if(prob) delete[] prob_;
   if(mean) delete[] mean_;
   if(stddev) delete[] stddev_;
   if(utility) { for(int j = 0; j < instance->getNbFacilities()*instance->getNbDisrType(); ++j) delete[] utility_[j]; }
   if(utility) delete[] utility_;
   if(sum_utility) delete[] sum_utility_; 
}

void SAAProblem::discreteChoiceDistribution(int i){
   // Facilities capacity (do not consider dummy facility)
   int nb_protec_node = instance->getNbFacilityProtLevels();
   int nb_sce = saainstance->getNbScenarios(i);
   for(int f = 0; f < instance->getNbFacilities(); ++f){
      int n = instance->getNodeFacility(f);
      for(int s = 0; s < nb_sce; ++s){  
         GRBLinExpr flow;
         for(int a: instance->getArcsArriving(n)){
            flow -= z[i][a*nb_sce+s];
         }
         for(int a: instance->getArcsLeaving(n)){
            flow += z[i][a*nb_sce+s];
         }
         for (int p = 0; p < nb_protec_node; ++p){
            flow -= x[i][f*nb_protec_node+p]*saainstance->getCapaFacilityProtection(i,s,f,p);
         }
         model[i]->addConstr(flow <= 0, "Flow_Facility_" + std::to_string(f) + 
                                        "_Scenario_" + std::to_string(s));
      }     
   }
}

void SAAProblem::normalDistribution(int i){
   int nb_disr_type = instance->getNbDisrType();
   int nb_protec_node = instance->getNbFacilityProtLevels();
   int nb_facilities = instance->getNbFacilities();
   int nb_sce = saainstance->getNbScenarios(i);
   
   mean[i] = model[i]->addVars(nb_facilities*nb_disr_type,GRB_CONTINUOUS);
   stddev[i] = model[i]->addVars(nb_facilities*nb_disr_type,GRB_CONTINUOUS);

   // Limite values for mean and std and naming variables
   for(int dt = 0; dt < nb_disr_type; ++dt) {
      for(int n = 0; n < nb_facilities; ++n){
         mean[i][dt*nb_facilities+n].set(GRB_DoubleAttr_LB,0.0);
         mean[i][dt*nb_facilities+n].set(GRB_StringAttr_VarName,
                                 "mean[" + std::to_string(dt) + 
                                 "," + std::to_string(n) + "]");

         stddev[i][dt*nb_facilities+n].set(GRB_DoubleAttr_LB,0.0);
         stddev[i][dt*nb_facilities+n].set(GRB_StringAttr_VarName,
                                 "stddev[" + std::to_string(dt) + 
                                 "," + std::to_string(n) + "]");
      }
   }

   // Mean should depend on the protection decisions
   double rho = 1 + (instance->getNbFacilities()-1)*instance->getFacilitiesImpact();
   for(int dt = 0; dt < nb_disr_type; ++dt) {
      for(int n = 0; n < nb_facilities; ++n){
         GRBLinExpr def_mean = rho*mean[i][dt*nb_facilities+n];
         GRBLinExpr def_stddev = rho*stddev[i][dt*nb_facilities+n];
         for(int nn = 0; nn < nb_facilities; ++nn){
            double weight = (nn == n) ? 1.0 : instance->getFacilitiesImpact();
            def_stddev -= (weight*instance->getStdDevNormalDist(dt,nn));
            for(int p = 0; p < nb_protec_node; ++p){
               def_mean -= (weight*instance->getMeanNormalDist(dt,nn,p)*x[i][nn*nb_protec_node+p]);
            }
         }
         model[i]->addConstr(def_mean == 0.0 , "Mean_Definition_Disr_" + std::to_string(dt) + 
                                               "_Facility/Node_" + std::to_string(n));
         model[i]->addConstr(def_stddev == 0.0 , "Stddev_Definition_Disr_" + std::to_string(dt) + 
                                               "_Facility/Node_" + std::to_string(n));
      }
   }
   
   // Facilities capacity (do not consider dummy facility)
   for(int f = 0; f < nb_facilities; ++f){
      int n = instance->getNodeFacility(f);
      for(int s = 0; s < nb_sce; ++s){  
         GRBLinExpr flow;
         for(int a: instance->getArcsArriving(n)){
            flow -= z[i][a*nb_sce+s];
         }
         for(int a: instance->getArcsLeaving(n)){
            flow += z[i][a*nb_sce+s];
         }
         int dt = saainstance->getScenarioDisruption(i,s);
         if(dt == -1){
            model[i]->addConstr(flow <= instance->getFullCapacity(), 
                                 "Flow_Facility_" + std::to_string(f) + 
                                 "_Scenario_" + std::to_string(s));
         }else{
            flow -= (mean[i][dt*nb_facilities+f] + stddev[i][dt*nb_facilities+f]*saainstance->getNormalDistRV(i,s,f));
            model[i]->addConstr(flow <= 0, "Flow_Facility_" + std::to_string(f) + 
                                           "_Scenario_" + std::to_string(s));
         }
      }     
   }
}

void SAAProblem::binomialDistribution(int i){
   int nb_scenarios = saainstance->getNbScenarios(i);
   int nb_disr_type = instance->getNbDisrType();
   int nb_capa_levels = instance->getNbBinomialCapaStates();
   int nb_facilities = instance->getNbFacilities();
   int nb_protec_node = instance->getNbFacilityProtLevels();

   prob[i] = model[i]->addVars(nb_facilities*nb_disr_type,GRB_CONTINUOUS);
   for(int dt = 0; dt < nb_disr_type; ++dt){
      for(int n = 0; n < nb_facilities; ++n){
         prob[i][dt*nb_facilities+n].set(GRB_DoubleAttr_LB,0.0);
         prob[i][dt*nb_facilities+n].set(GRB_DoubleAttr_UB,1.0);
      }
   }

   capa[i] = new GRBVar*[nb_scenarios];
   for(int s = 0; s < nb_scenarios; ++s){
      capa[i][s] = model[i]->addVars(nb_facilities*nb_capa_levels,GRB_BINARY);
   }

   // Capacities definition
   for(int n = 0; n < nb_facilities; ++n){
      for(int w = 0; w < nb_capa_levels; ++w){
         for(int s = 0; s < nb_scenarios; ++s){
            int dt = saainstance->getScenarioDisruption(i,s);
            if(dt != -1){
               GRBLinExpr exp1 = saainstance->getBinomialUniformRV(i,s,n,w)*capa[i][s][n*nb_capa_levels+w]
                                 - prob[i][dt*nb_facilities+n] + instance->getEpsilon();
               GRBLinExpr exp2 = capa[i][s][n*nb_capa_levels+w] - prob[i][dt*nb_facilities+n] + saainstance->getBinomialUniformRV(i,s,n,w);
               model[i]->addConstr(exp1 <= 0, "1_Capacity_Facility/Node_" + std::to_string(n) + 
                                             "_Level_" + std::to_string(w) + 
                                             "_Scenario" + std::to_string(s));
               model[i]->addConstr(exp2 >= 0, "2_Capacity_Facility/Node_" + std::to_string(n) + 
                                             "_Level_" + std::to_string(w) + 
                                             "_Scenario" + std::to_string(s));
               
               //Trying with indicator constraint
               model[i]->addGenConstrIndicator(capa[i][s][n*nb_capa_levels+w], 1, 
                              prob[i][dt*nb_facilities+n] >= saainstance->getBinomialUniformRV(i,s,n,w) + instance->getEpsilon(),
                              "1Implication_Scenario_" + std::to_string(s) +
                              "Facility/node" + std::to_string(n) +
                              "_Level_" + std::to_string(w));
               model[i]->addGenConstrIndicator(capa[i][s][n*nb_capa_levels+w], 0, 
                              prob[i][dt*nb_facilities+n] <= saainstance->getBinomialUniformRV(i,s,n,w),
                              "2Implication_Scenario_" + std::to_string(s) +
                              "Facility/node" + std::to_string(n) +
                              "_Level_" + std::to_string(w));
            }
         }
      }
   }

   // Ordered constraint
   if(useinequality){
      for(int dt = 0; dt < nb_disr_type; ++dt){
         for(int f = 0; f < nb_facilities; ++f){
            std::vector<std::tuple<int,int>> ordered_index = saainstance->getOrderedBernoulliUniformRV(i,dt,f);
            if(ordered_index.size() >= 2){
               for(int j = 0; j < (int)ordered_index.size()-1; ++j){
                  int s = std::get<0>(ordered_index[j]);
                  int w = std::get<1>(ordered_index[j]);
                  int ss = std::get<0>(ordered_index[j+1]);
                  int ww = std::get<1>(ordered_index[j+1]);
                  model[i]->addConstr(capa[i][s][f*nb_capa_levels+w] >= 
                                    capa[i][ss][f*nb_capa_levels+ww],
                                    "Ordered_Capacity_" + std::to_string(j) +
                                    "_Disruption_" + std::to_string(dt) +
                                    "_Node/facility_" + std::to_string(f));
               }
            }
         }
      }
   }

   double rho = 1 + (nb_facilities-1)*instance->getFacilitiesImpact();
   // Probability equal to weighted average
   for(int dt = 0; dt < nb_disr_type; ++dt){
      for(int n = 0; n < nb_facilities; ++n){
         GRBLinExpr expprob = rho*prob[i][dt*nb_facilities+n];
         for(int nn = 0; nn < nb_facilities; ++nn){
            double weight = (nn == n) ? 1.0 : instance->getFacilitiesImpact();
            for(int p = 0; p < nb_protec_node; ++p){
               expprob -= weight*instance->getBernoulliProbCapa(dt,nn,p)*x[i][nn*nb_protec_node+p];
            }
         }
         model[i]->addConstr(expprob == 0, "Probability_Definition_Disruption_" + std::to_string(dt) + "_Facility_" + std::to_string(n));
      }
   }

   // Facilities capacity (do not consider dummy facility)
   for(int f = 0; f < nb_facilities; ++f){
      int n = instance->getNodeFacility(f);
      for(int s = 0; s < nb_scenarios; ++s){  
         GRBLinExpr flow;
         for(int a: instance->getArcsArriving(n)){
            flow -= z[i][a*nb_scenarios+s];
         }
         for(int a: instance->getArcsLeaving(n)){
            flow += z[i][a*nb_scenarios+s];
         }
         int dt = saainstance->getScenarioDisruption(i,s);
         if(dt == -1){
            model[i]->addConstr(flow <= instance->getFullCapacity(), 
                                 "Flow_Facility_" + std::to_string(f) + 
                                 "_Scenario_" + std::to_string(s));
         }else{
            for(int w = 0; w < nb_capa_levels; ++w){
               flow -= capa[i][s][f*nb_capa_levels+w]*instance->getCapaPerStage();
            }
            model[i]->addConstr(flow <= 0, "Flow_Facility_" + std::to_string(f) + 
                                           "_Scenario_" + std::to_string(s));
         }
      }     
   }
}

void SAAProblem::stdNormaDistribution(int i){
   int nb_scenarios = saainstance->getNbScenarios(i);
   int nb_disr_type = instance->getNbDisrType();
   int nb_capa_levels = instance->getNbCapaStates();
   int nb_protec_node = instance->getNbFacilityProtLevels();
   int nb_facilities = instance->getNbFacilities();

   // Defining variables
   utility[i] = new GRBVar*[nb_facilities*nb_disr_type];
   sum_utility[i] = model[i]->addVars(nb_facilities*nb_disr_type,GRB_CONTINUOUS);
   for(int dt = 0; dt < nb_disr_type; ++dt){
      for(int n = 0; n < nb_facilities; ++n){
         utility[i][dt*nb_facilities+n] = model[i]->addVars(nb_capa_levels,GRB_CONTINUOUS);
         for(int w = 0; w < nb_capa_levels; ++w){
            utility[i][dt*nb_facilities+n][w].set(GRB_DoubleAttr_LB,0.0);
            utility[i][dt*nb_facilities+n][w].set(GRB_StringAttr_VarName,
                                 "u[" + std::to_string(dt) + 
                                 "," + std::to_string(n) + 
                                 "," + std::to_string(w) + "]");
         }
         sum_utility[i][dt*nb_facilities+n].set(GRB_StringAttr_VarName,
                                 "sum_u[" + std::to_string(dt) + 
                                 "," + std::to_string(n) + "]");
      }
   }
   
   capa[i] = new GRBVar*[nb_scenarios];
   for(int s = 0; s < nb_scenarios; ++s){
      capa[i][s] = model[i]->addVars(nb_facilities*nb_capa_levels,GRB_BINARY);
      for(int n = 0; n < nb_facilities; ++n){
         for(int w = 0; w < nb_capa_levels; ++w){
            capa[i][s][n*nb_capa_levels+w].set(GRB_StringAttr_VarName,
                                 "capa[" + std::to_string(s) + 
                                 "," + std::to_string(n) + 
                                 "," + std::to_string(w) + "]");
         }
      }
   }

   cum_capa[i] = new GRBVar*[nb_scenarios];
   for(int s = 0; s < nb_scenarios; ++s){
      cum_capa[i][s] = model[i]->addVars(nb_facilities*nb_capa_levels,GRB_BINARY);
      for(int n = 0; n < nb_facilities; ++n){
         for(int w = 0; w < nb_capa_levels; ++w){
            cum_capa[i][s][n*nb_capa_levels+w].set(GRB_StringAttr_VarName,
                                 "cum_capa[" + std::to_string(s) + 
                                 "," + std::to_string(n) + 
                                 "," + std::to_string(w) + "]");
         }
      }
   }

   // Capacity definition
   for(int s = 0; s < nb_scenarios; ++s){
      for(int n = 0; n < nb_facilities; ++n){
         // Choose exactly one realization/capacity
         GRBLinExpr sum_capa;
         for(int w = 0; w < nb_capa_levels; ++w){
            sum_capa += capa[i][s][n*nb_capa_levels+w];
         }
         model[i]->addConstr(sum_capa == 1, "Sum_Capacity_Scenario_" + std::to_string(s) + 
                                            "_Node/Facility_"+ std::to_string(n));

         // Define partial sum of capacities
         model[i]->addConstr(cum_capa[i][s][n*nb_capa_levels+(nb_capa_levels-1)] == 1, 
                                                   "Auxiliar_Capacity_+" + std::to_string(nb_capa_levels-1) + 
                                                   "_Scenario_" + std::to_string(s) + 
                                                   "_Node/Facility_"+ std::to_string(n));

         for(int w = 0; w < (nb_capa_levels-1); ++w){
            GRBLinExpr cumsum_capa = cum_capa[i][s][n*nb_capa_levels+w] - cum_capa[i][s][n*nb_capa_levels+(w+1)] + capa[i][s][n*nb_capa_levels+(w+1)];
            model[i]->addConstr(cumsum_capa == 0, "Auxiliar_Capacity_+" + std::to_string(w) + 
                                                  "_Scenario_" + std::to_string(s) + 
                                                  "_Node/Facility_"+ std::to_string(n));
         }
      }
   }
   
   // Utility definition
   for(int dt = 0; dt < nb_disr_type; ++dt){
      for(int n = 0; n < nb_facilities; ++n){
         GRBLinExpr exp_sum_utility = sum_utility[i][dt*nb_facilities+n];
         for(int w = 0; w < nb_capa_levels; ++w){
            exp_sum_utility -= utility[i][dt*nb_facilities+n][w];
            
            GRBLinExpr exp_utility = utility[i][dt*nb_facilities+n][w];
            for(int p = 0; p < nb_protec_node; ++p){
               exp_utility -= instance->getUtility(dt,n,p,w)*x[i][n*nb_protec_node+p];
            }
            for(int n2 = 0; n2 < nb_facilities; ++n2){
               if(n2 != n){
                  for(int p = 0; p < nb_protec_node; ++p){
                     exp_utility -= instance->getFacilitiesImpact()*
                                    instance->getUtility(dt,n2,p,w)*x[i][n2*nb_protec_node+p];
                  }
               }
            }
            model[i]->addConstr(exp_utility == 0,"Utility_Definition_Disrup_" + std::to_string(dt) + 
                                                 "_Facility/Node_" + std::to_string(n) +
                                                 "_Capacity_" + std::to_string(w));
         } 
         model[i]->addConstr(exp_sum_utility == 0,"Sum_Utility_Definition_Disrup_" + std::to_string(dt) + 
                                                 "_Facility/Node_" + std::to_string(n));
      }
   }

   // McCormick
   mccormick1[i] = new GRBVar*[nb_scenarios];
   for(int s = 0; s < nb_scenarios; ++s){
      mccormick1[i][s] = model[i]->addVars(nb_facilities*(nb_capa_levels-1),GRB_CONTINUOUS);
      for(int n = 0; n < nb_facilities; ++n){
         for(int w = 0; w < nb_capa_levels-1; ++w){
            mccormick1[i][s][n*(nb_capa_levels-1)+w].set(GRB_StringAttr_VarName,
                                 "mccormick[" + std::to_string(s) + "," 
                                 + std::to_string(n) + 
                                 "," + std::to_string(w) + "]");
         }
      }
   }

   for(int s = 0; s < nb_scenarios; ++s){
      int dt = saainstance->getScenarioDisruption(i,s);
      if(dt != -1){
         for(int n = 0; n < nb_facilities; ++n){
            for(int w = 0; w < (nb_capa_levels-1); ++w){
               GRBLinExpr exp1 = mccormick1[i][s][n*(nb_capa_levels-1)+w] - sum_utility[i][dt*nb_facilities+n];
               GRBLinExpr exp2 = mccormick1[i][s][n*(nb_capa_levels-1)+w] - instance->getMaxUtility(dt,n)*cum_capa[i][s][n*nb_capa_levels+w];
               GRBLinExpr exp3 = mccormick1[i][s][n*(nb_capa_levels-1)+w] - sum_utility[i][dt*nb_facilities+n] 
                              + instance->getMaxUtility(dt,n) - instance->getMaxUtility(dt,n)*cum_capa[i][s][n*nb_capa_levels+w];

               model[i]->addConstr(exp1 <= 0,"McCormick1_Scenario_"+ std::to_string(s) +
                                             "_Disrup_" + std::to_string(dt) + 
                                             "_Facility/Node_" + std::to_string(n) +
                                             "_Capacity_" + std::to_string(w));
               model[i]->addConstr(exp2 <= 0,"McCormick2"+ std::to_string(s) +
                                             "_Disrup_"+  std::to_string(dt) + 
                                             "_Facility/Node_" + std::to_string(n) +
                                             "_Capacity_" + std::to_string(w));
               model[i]->addConstr(exp3 >= 0,"McCormick3"+ std::to_string(s) +
                                             "_Disrup_"+ std::to_string(dt) + 
                                             "_Facility/Node_" + std::to_string(n) +
                                             "_Capacity_" + std::to_string(w));
            }
         }
      }
   }

   // Final capa definition
   for(int s = 0; s < nb_scenarios; ++s){
      int dt = saainstance->getScenarioDisruption(i,s);
      if(dt != -1){
         for(int n = 0; n < nb_facilities; ++n){
            for(int w = 0; w < (nb_capa_levels-1); ++w){
               GRBLinExpr capa = saainstance->getStdNormaUniformRV(i,s,n)*mccormick1[i][s][n*(nb_capa_levels-1)+w];
               for(int j = 0; j <= w; ++j){
                  capa -= utility[i][dt*nb_facilities+n][j];
               }
               model[i]->addConstr(capa <= (-instance->getEpsilon()), "Def1_Capacity_" + std::to_string(w) + "_Scenario_"+ std::to_string(s) + "_Node/Facility_"+ std::to_string(n));
            }
            for(int w = 1; w < nb_capa_levels; ++w){
               GRBLinExpr capa = (1.0-saainstance->getStdNormaUniformRV(i,s,n))*sum_utility[i][dt*nb_facilities+n];
               capa -= (1.0-saainstance->getStdNormaUniformRV(i,s,n))*mccormick1[i][s][n*(nb_capa_levels-1)+(w-1)];
               for(int j = w; j < nb_capa_levels; ++j){
                  capa -= utility[i][dt*nb_facilities+n][j];
               }
               model[i]->addConstr(capa <= 0, "Def2_Capacity_" + std::to_string(w) + "_Scenario_"+ std::to_string(s) + "_Node/Facility_"+ std::to_string(n));
            }           
         }
      }
   }

   // Ordered constraint
   if(useinequality){
      for(int dt = 0; dt < nb_disr_type; ++dt){
         for(int n = 0; n < nb_facilities; ++n){
            std::vector<int> ordered_index = saainstance->getOrderedStdNormaUniformRV(i,dt,n);
            if(ordered_index.size() >= 2){
               for(int j = 0; j < (int)ordered_index.size()-1; ++j){
                  int s = ordered_index[j];
                  int ss = ordered_index[j+1];
                  for(int w = 0; w < nb_capa_levels; ++w){
                     GRBLinExpr ordered_capa = cum_capa[i][s][n*nb_capa_levels+w] 
                                             - cum_capa[i][ss][n*nb_capa_levels+w];
                     model[i]->addConstr(ordered_capa >= 0,
                                       "Ordered_" + std::to_string(j) +
                                       "_Disruption_" + std::to_string(dt) +
                                       "_Node/facility_" + std::to_string(n) +
                                       "_Capacity_" + std::to_string(w));
                  }
               }
            }
         }
      }
   }

   // Facilities capacity (do not consider dummy facility)
   for(int f = 0; f < nb_facilities; ++f){
      int n = instance->getNodeFacility(f);
      for(int s = 0; s < nb_scenarios; ++s){  
         GRBLinExpr flow;
         for(int a: instance->getArcsArriving(n)){
            flow -= z[i][a*nb_scenarios+s];
         }
         for(int a: instance->getArcsLeaving(n)){
            flow += z[i][a*nb_scenarios+s];
         }
         int dt = saainstance->getScenarioDisruption(i,s);
         if(dt == -1){
            model[i]->addConstr(flow <= instance->getFullCapacity(), 
                                 "Flow_Facility_" + std::to_string(f) + 
                                 "_Scenario_" + std::to_string(s));
         }else{
            for(int w = 0; w < nb_capa_levels; ++w){
               flow -= capa[i][s][f*nb_capa_levels+w]*instance->getCapacity(w);
            }
            model[i]->addConstr(flow <= 0, "Flow_Facility_" + std::to_string(f) + 
                                           "_Scenario_" + std::to_string(s));
         }
      }     
   }
}

void SAAProblem::generateCommonVariables(int i){
   int nb_scenarios = saainstance->getNbScenarios(i);
   int nb_protec_node = instance->getNbFacilityProtLevels();

   count_x = instance->getNbFacilities()*nb_protec_node;
   count_y = instance->getNbEdges();
   int count_z = instance->getNbArcsWithDummy()*nb_scenarios;

   x[i] = model[i]->addVars(count_x,GRB_BINARY);
   y[i] = model[i]->addVars(count_y,GRB_BINARY);
   z[i] = model[i]->addVars(count_z,GRB_CONTINUOUS);

   if(instance->getConstOrObjCost() == true){
      for (int f = 0; f < instance->getNbFacilities(); ++f){
         for(int p = 0; p < nb_protec_node; ++p){
            x[i][f*nb_protec_node+p].set(GRB_DoubleAttr_Obj,instance->getCostFacilityProtection(f,p));
         }
      }
      for (int e = 0; e < instance->getNbEdges(); ++e){
         y[i][e].set(GRB_DoubleAttr_Obj,instance->getCostEdge(e));
      }
   } 

   // Naming variables
   for (int f = 0; f < instance->getNbFacilities(); ++f){
      for(int p = 0; p < nb_protec_node; ++p) {
         x[i][f*nb_protec_node+p].set(GRB_StringAttr_VarName,"x[" + std::to_string(f) + "," + std::to_string(p) + "]");
      }
   }
   for (int e = 0; e < instance->getNbEdges(); ++e) y[i][e].set(GRB_StringAttr_VarName,"y[" + std::to_string(e) + "]");
   for(int a = 0; a < instance->getNbArcsWithDummy(); ++a){
      double length = instance->getArcLength(a);
      std::pair<int,int> arc = instance->getArc(a);
      for(int s = 0; s < nb_scenarios; ++s){
         z[i][a*nb_scenarios+s].set(GRB_DoubleAttr_Obj,length*saainstance->getWeightScenarioSAA(i,s));
         z[i][a*nb_scenarios+s].set(GRB_StringAttr_VarName,"z[" + std::to_string(s) + 
                                                      "," + std::to_string(arc.first) + 
                                                      "," + std::to_string(arc.second) + "]");
      }
   }
   for(int a = 0; a < instance->getNbArcs(); ++a){
      for(int s = 0; s < nb_scenarios; ++s) 
         z[i][a*nb_scenarios+s].set(GRB_DoubleAttr_UB,instance->getArcMaximumCapacity());
   }
}

void SAAProblem::generateCommonConstraints(int i){
   int nb_faclities = instance->getNbFacilities();
   int nb_protec_node = instance->getNbFacilityProtLevels();

   // Facilities protection
   for (int f = 0; f < nb_faclities; ++f){
      GRBLinExpr invest = 0;
      for (int p = 0; p < nb_protec_node; ++p){
         invest += x[i][f*nb_protec_node+p];
      }
      model[i]->addConstr(invest == 1, "Investment_Facility_" + std::to_string(f));
   }

   // Budget
   if(instance->getConstOrObjCost() == false){
      GRBLinExpr cost = 0;
      for (int f = 0; f < nb_faclities; ++f){
         for(int p = 0; p < nb_protec_node; ++p){
            cost += x[i][f*nb_protec_node+p]*instance->getCostFacilityProtection(f,p);
         }
      }
      for (int e = 0; e < instance->getNbEdges(); ++e){
         cost += y[i][e]*instance->getCostEdge(e);
      }
      model[i]->addConstr(cost <= instance->getBudget(), "Budget");
   }

   // Constraints for second stage. 
   int nb_scenarios = saainstance->getNbScenarios(i);

   // Clients demand 
   for(int c = 0; c < instance->getNbClients(); ++c){
      double demand = instance->getClientDemand(c);
      int n = instance->getNodeClient(c);
      for(int s = 0; s < nb_scenarios; ++s){
         GRBLinExpr flow;
         for(int a: instance->getArcsArriving(n)){
            flow += z[i][a*nb_scenarios+s];
         }
         for(int a: instance->getArcsLeaving(n)){
            flow -= z[i][a*nb_scenarios+s];
         }
         model[i]->addConstr(flow == demand, "Flow_Client_" + std::to_string(c) + 
                                             "_Scenario_" + std::to_string(s));
      }  
   }

   // Arcs capacity (do not consider dummy facilities)
   for(int a = 0; a < instance->getNbArcs(); ++a){
      for(int s = 0; s < nb_scenarios; ++s){
         int e = instance->getEdgeFromArc(a);
         GRBLinExpr capacity = z[i][a*nb_scenarios+s] - instance->getArcMaximumCapacity()*y[i][e];
         model[i]->addConstr(capacity <= 0, "Capacity_Arc_" + std::to_string(a) + 
                                            "_Scenario_" + std::to_string(s));
      }
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////