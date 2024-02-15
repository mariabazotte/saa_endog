#ifndef INSTANCE_HPP
#define INSTANCE_HPP

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <math.h>
#include <random>
#include "../input/input.hpp"


class Instance {
    protected:
        const Input & input;                   /* Input info. */
        std::string network_file;              /* File with network information. */
        const std::string parameters_file;     /* File with parameters information. */
        const double seed;                     /* Seed for data generation. */
        std::string name_instance;             /* Name instance. */

        bool UNCERTAIN_NO_DISR = true;         /* True if no disruption is also uncertain. */
        bool ALL_DISR_LEVELS   = true;         /* True if facilities in the disruption level 0 also are impacted by disruption (always true).*/

        int nb_nodes;           /* Number of nodes. */
        int nb_edges;           /* Number of edges. */
        int nb_arcs;            /* Number of arcs. */
        int nb_facilities;      /* Number of installed facilities.*/
        int nb_clients;         /* Number of client nodes. */
        int nb_disr_type;       /* Number of disruption types. e.g earthquake, snow storm */
        int nb_disr_level;      /* Number of levels/intensity of disruptions. Usually is 3 levels (0,1,2) or (1,2,3). 0 represents no impact. */

        int nb_prot_level_facility;        /* Number of protection levels for facilities. */
        int nb_capa_states;                /* Number of possible capacity states/values after disruption. */
        double penalty_demand;             /* Penalty for each demand that is not served. */
        bool const_or_obj_cost;            /* false: constraint with limit budget for facility and edge installation, true: minimize the installation cost. */
        double percentage_budget;          /* Percentage to compute the maximum budget. */
        double budget;                     /* Maximum budget for installation cost when const_or_obj_cost = false. */
        double facilities_impact;          /* Impact factor of facilities on the failure of other facilities. */
        double max_flow_arcs;              /* Maximum flow of arcs/edges. */
        int num_dummy_facility;            /* Id node for dummy facility. */
        int num_rep_edge_dummy_facility;   /* Number of representative edge for dummy facility. */
        int nb_arcs_dummy;                 /* Number arcs related to dummy facility. */

        std::map<std::string, int> id_cities;           /* Node id corresponding to the city. */
        std::vector<std::string> name_cities;           /* City corresponding to node id. */
        std::vector<int> node_demand;                   /* Demand of nodes. (general formulation) */
        std::vector<int> client_demand;                 /* Demand of clients. (simple formulation) */
        std::vector<double> edge_length;                /* Length of edges. */
        std::vector<std::vector<int>> node_disr_level;  /* node_disr_level[node][disr_type] : level of disruption at node considering specific disruption type. */
        std::vector<int> facilities;                    /* Node id's of installed facilities. */
        std::vector<int> clients;                       /* Node id's of client nodes. */

        std::vector<std::vector<int>> adj_edges;           /* adjacent_edges[n] : list of edges adjacent to node n. */
        std::vector<int> edge_of_arc;                      /* Edge that arc represents. */
        std::vector<std::pair<int, int>> arcs;             /* <source node, target node> of arcs. */
        std::vector<std::vector<int>> arcs_leaving_node;   /* arcs_leaving_node[n] : list of arcs leaving node n. */
        std::vector<std::vector<int>> arcs_arriving_node;  /* arcs_leaving_node[n] : list of arcs arriving at node n. */    

        std::vector<std::string> cities_disr;  /* Local/center each disruption type occurs. */
        std::vector<double> prob_disr_type;    /* Probability of each disruption type to happen.*/
        double prob_no_disr;                   /* 1-sum(prob_disr_type) = probability normal case. */

        int max_capacity;                                     /* Value of maximum capacity for facilities. */
        double capa_per_state;                                /* Increase of capacity according to capacity level/state. */
        double max_multi_inst_edge;                           /* Value of maximum multiplicator for installing edge cost. */
        std::vector<double> capacity;                         /* capacity[state] : capacity value of nodes acording to possible states. */
        double max_capacity_normal;                           /* Maximum capacity for normal distribution (continuous variable). */
        std::vector<std::vector<double>> cost_inst_facility;  /* cost_inst_facility[n][p] : cost installing facility at node n with protection p.*/
        std::vector<double> max_cost_inst_facility;           /* max_cost_inst_facility[n]: Max cost install facility node n.*/

        std::vector<std::vector<std::vector<double>>> bernoulli_prob_capa;              /* bernoulli_prob_capa[t][n][p] : Probability of success considering dirspution type t, node n and protection p. */
        std::vector<std::vector<std::vector<std::vector<double>>>> binomial_prob_capa;  /* binomial_prob_capa[t][n][p][l] : Probability of capacity level l considering disruption type t, node n and protection level p. */
        std::vector<std::vector<std::vector<double>>> bernoulli_prob_not_fail;          /* bernoulli_prob_not_fail[t][e][p] : Probability of edge e not failing in disruption type t and with protection level p. (prob of failing : 1-prob_fail[t][e][p]). */

        std::vector<std::vector<std::vector<std::vector<double>>>> utility;             /* utility[t][n][p][l]: utility for disruption t, node n, protection p and capacity level l.*/
        std::vector<std::vector<double>> max_utility;                                   /* max_utility[t][n]: max utility for disruption t and node n for standard normalization*/
        std::vector<std::vector<double>> mean_capa;                                     /* mean_capa[node][protec]: mean capacity for node n and protection p. */

        std::vector<std::vector<std::vector<double>>> mean_normal_dist;                 /* mean_normal_dist[dt][n][p]: mean parameter normal distribution for disruption dt, node n and prtection p. */
        std::vector<std::vector<double>> stddev_normal_dist;                            /* stddev_normal_dist[dt][n]: std deviation parameter normal distribution for disruption dt, node n. */

        std::string getValueFromParameterFile(std::ifstream &, std::string);            /* Function to read each value from parameters file. */
        void readNetworkFile();                                                         /* Function to read network file.*/
        void readParametersFile();                                                      /* Function to read parameters file.*/
        void generateData();                                                            /* Dunction to generate data. */
        void computeProbabilities();
        void computeUtilities();
        void computeFacilitiesProbabilities();
        void computeFacilitiesMean();
        void computeParamsNormalDistr();
    
    public:
        Instance(const Input & input);
        ~Instance() {}

        std::string getNetworkFile() const { return network_file; }
        std::string getParametersFile() const { return parameters_file; }
        std::string getInstanceName() const { return name_instance; }
        std::string getNodeName(int n) const { return name_cities[n]; }
        std::string getArcName(int a) { return (name_cities[arcs[a].first] + "," + name_cities[arcs[a].second]); }

        double getEpsilon() const { return 0.0000001;}

        int getNbNodes() const { return nb_nodes; }
        int getNbEdges() const { return nb_edges; }
        int getNbArcs() const { return nb_arcs; }
        int getNbNodesWithDummy() const { return nb_nodes + 1; }
        int getNbArcsWithDummy() const { return nb_arcs + nb_arcs_dummy; }
        int getDummyFacility() const { return num_dummy_facility; }
        int getNbDisrType() const { return nb_disr_type; }
        int getNbDisrLevel() const { return nb_disr_level; }
        int getNbCapaStates() const { return nb_capa_states+1; }
        int getNbBinomialCapaStates() const { return nb_capa_states; }
        int getNbFacilityProtLevels() const { return nb_prot_level_facility; }
        int getNbFacilities() const { return nb_facilities; }

        double getCapaPerStage() const { return capa_per_state; }
        double getCapacity(int l) const { return capacity[l]; }
        double getFullCapacity() const { return (double)max_capacity; }
        double getFacilitiesImpact() const { return facilities_impact; }
        double getBudget() const { return budget; }
        bool getConstOrObjCost() const { return const_or_obj_cost; }
        
        const std::pair<int, int> & getArc(int a) const { return arcs[a]; }
        int getEdgeFromArc(int a) const { return edge_of_arc[a]; }
        double getArcLength(int a) const { return edge_length[edge_of_arc[a]]/100; }
        double getArcMaximumCapacity() const { return max_flow_arcs; }

        const std::vector<int> & getArcsLeaving(int n) const { return arcs_leaving_node[n]; }
        const std::vector<int> & getArcsArriving(int n) const { return arcs_arriving_node[n]; }

        const std::vector<double> & getProbDisrType() const { return prob_disr_type; } 
        double getProbDisrType(int dt) const { return prob_disr_type[dt]; }
        double getProbNoDisr() const { return prob_no_disr; }
        bool getUncertainNoDisr() const { return UNCERTAIN_NO_DISR; }
        double getEdgeCostMultiplier() const { return max_multi_inst_edge; }

        int getNbFacilitiesWithDummy() const { return nb_facilities + 1; }
        int getNbClients() const { return nb_clients; }
        int getNodeClient(int c) const { return clients[c]; }
        int getNodeFacility(int f) const { return facilities[f]; }
        int getClientIndex(int n) const { return std::distance(clients.begin(),std::find(clients.begin(),clients.end(),n)); }
        int getFacilityIndex(int n) const { return std::distance(facilities.begin(),std::find(facilities.begin(),facilities.end(),n)); }

        int getClientDemand(int c) const { return client_demand[c]; }
        bool getIsFacility(int n) const  { return std::find(facilities.begin(), facilities.end(), n) != facilities.end() ;}
        double getCostFacilityProtection(int f, int p) const { return cost_inst_facility[f][p]; }
        double getCostEdge(int e) const { return max_multi_inst_edge*edge_length[e]; }

        double getMaxUtility(int dt, int f) const { return max_utility[dt][f]; }
        double getUtility(int dt, int f,int p, int w) const { return utility[dt][f][p][w]; }
        double getMeanCapacity(int f,int p) const { return mean_capa[f][p]; }
        double getMeanNormalDist(int dt, int f, int p) const { return mean_normal_dist[dt][f][p]; }
        double getStdDevNormalDist(int dt, int f) const {return stddev_normal_dist[dt][f]; }
        double getBernoulliProbCapa(int dt,int f,int p) const { return bernoulli_prob_capa[dt][f][p]; }
        const std::vector<double> & getBinomialProbCapaVector(int dt,int f,int p) const { return binomial_prob_capa[dt][f][p]; }
        const std::vector<int> & getAdjacentEdges(int f) const { return adj_edges[facilities[f]]; } 
        double getRhsCut(const double* sol) const { return std::inner_product(client_demand.begin(), client_demand.end(), sol, 0.0); } 

        std::string write() const;
};

#endif
