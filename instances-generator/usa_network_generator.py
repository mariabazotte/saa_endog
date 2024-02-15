import utm
import pyproj
from scipy.spatial import Delaunay
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import math
import argparse
import gurobipy as gp
from gurobipy import GRB

# Alabama, Arkansas, Florida, Georgia, 
# Kentucky, Louisiana, Mississippi, 
# North Carolina, South Carolina, Tennessee
states = ['AL','AR','FL','GA','KY','LA','MS','NC','SC','TN']

# 3 different disruption types
cities_disr = ['Tampa','Raleigh','Huntsville']
prob_disr = [0.1,0.05,0.1]
nb_disr_type = 3

# 3 different levels of disruption (radius distance to decide disruption level for each disruption type)
distance_levels = [[250,800],[160,800],[80,400]]
nb_disr_level = 3

multiplier_transp_cost = 10
median_house_cost = 339048 #zillow 2023

# 2000000 -> 3 nodes
# 1500000 -> 5 nodes
# 1000000 -> 10 nodes
# 650000 -> 15 nodes
# 500000 -> 19 nodes
# 403500 -> 25 nodes
# 350000 -> 30 nodes
# 300000 -> 39 nodes
# 250000 -> 48 nodes
# 200000 -> 48 nodes

def parse_arg():
    parser = argparse.ArgumentParser()
    parser.add_argument('-minpop', metavar='[minpop]', \
                        help = 'minimum population to select cities',\
                        type = int, default=1000000)
    parser.add_argument('-nbfacilities', metavar='[nbfacilities]', \
                        help = 'Number of facilities',\
                        type = int, default=3)
    arg = parser.parse_args()
    return arg

class Data():
    def __init__(self,arg):
        self.readCitiesFile(arg.minpop)
        self.computeDelaunay()
        self.computeDisruption()
        self.solveModel(arg.nbfacilities)
        self.writeFile()
    
    def readCitiesFile(self,minpop):
        self.min_population = minpop
        df = pd.read_csv('data_usa/uscities.csv')
        df = df[['city','state_id','lat','lng','population']]
        df = df.loc[df['state_id'].isin(states)]
        self.cities = df.loc[df['population'] >= self.min_population] 
        self.num_cities = self.cities.shape[0]
        self.cities_and_disr_places = pd.concat([self.cities,df.loc[(df['city'] == 'Tampa') & (df['state_id'] == 'FL')]])
        self.cities_and_disr_places = pd.concat([self.cities_and_disr_places,df.loc[(df['city'] == 'Raleigh') & (df['state_id'] == 'NC')]])
        self.cities_and_disr_places = pd.concat([self.cities_and_disr_places,df.loc[(df['city'] == 'Huntsville') & (df['state_id'] == 'AL')]])
        print(self.cities_and_disr_places)

    def computeDelaunay(self):
        lat = np.array(self.cities_and_disr_places['lat'])
        lon = np.array(self.cities_and_disr_places['lng'])
        proj_wgs84 = pyproj.Proj(init="epsg:4326")
        proj_ = pyproj.Proj(init="epsg:9311")
        x, y = pyproj.transform(proj_wgs84, proj_, lon, lat)
        self.points = np.array([[x[i],y[i]] for i in range(self.num_cities)])
        self.points_disr = np.array([[x[self.num_cities+i],y[self.num_cities+i]] for i in range(nb_disr_type)])
        
        tri = Delaunay(self.points)
        #plt.triplot(points[:,0], points[:,1], tri.simplices)
        #plt.plot(points[:,0], points[:,1], 'o')
        #plt.show()

        self.edges = set()
        for n in range(tri.nsimplex):
            edge = sorted([tri.vertices[n,0], tri.vertices[n,1]])
            self.edges.add((edge[0], edge[1]))
            edge = sorted([tri.vertices[n,0], tri.vertices[n,2]])
            self.edges.add((edge[0], edge[1]))
            edge = sorted([tri.vertices[n,1], tri.vertices[n,2]])
            self.edges.add((edge[0], edge[1]))

        for a,b in self.edges:
            dist =  math.pow((self.points[a,0] - self.points[b,0])**2 + (self.points[a,1] - self.points[b,1])**2, 1/2)
            print('dist ', a, ' ', b, ' = ', dist)
            
        graph = nx.Graph(list(self.edges))
        nx.draw(graph)
        #print(graph.edges())
        # plot graph
        pointIDXY = dict(zip(range(len(self.points)), self.points))
        nx.draw(graph, pointIDXY)
        plt.show()
    
    def computeDisruption(self):  
        self.disr1 = list()
        self.disr2 = list()
        self.disr3 = list()
        for i in range(self.num_cities):
            # first disruption type
            dist =  math.pow((self.points[i,0] - self.points_disr[0,0])**2 + (self.points[i,1] - self.points_disr[0,1])**2, 1/2)/1000
            if dist <= distance_levels[0][0] :
                self.disr1.append(2)
            elif dist <= distance_levels[0][1]:
                self.disr1.append(1)
            else:
                self.disr1.append(0)
            # second disruption type 
            dist =  math.pow((self.points[i,0] - self.points_disr[1,0])**2 + (self.points[i,1] - self.points_disr[1,1])**2, 1/2)/1000
            if dist <= distance_levels[1][0] :
                self.disr2.append(2)
            elif dist <= distance_levels[1][1]:
                self.disr2.append(1)
            else:
                self.disr2.append(0)
            # third disruption type
            dist =  math.pow((self.points[i,0] - self.points_disr[2,0])**2 + (self.points[i,1] - self.points_disr[2,1])**2, 1/2)/1000
            if dist <= distance_levels[2][0] :
                self.disr3.append(2)
            elif dist <= distance_levels[2][1]:
                self.disr3.append(1)
            else:
                self.disr3.append(0)
    
    def solveModel(self,nbfacilities):  
        self.nb_nodes = self.num_cities
        self.nb_edges = len(self.edges)
        self.nb_facilities = min(self.num_cities,nbfacilities)
        
        demand = list()
        for i in range(self.nb_nodes):
            demand.append(round(self.cities['population'].iloc[i]/10000))
        max_node_capacity = round(sum(demand)/(self.nb_facilities*0.9))
        max_arc_capacity = round(sum(demand)) # no restriction for arcs, big M
        
        arcs = list()
        dist_arcs = list()
        for a,b in self.edges:
            arcs.append((a,b))
            arcs.append((b,a))
            dist =  math.pow((self.points[a,0] - self.points[b,0])**2 + (self.points[a,1] - self.points[b,1])**2, 1/2)/1000
            dist_arcs.append(dist)
            dist_arcs.append(dist)
        nb_arcs = len(arcs)
    
        # Solving model to define the facilities
        
        # Melkote, S., & Daskin, M. S. (2001). Capacitated facility location/network design problems. 
        # European journal of operational research, 129(3), 481-495.
        
        model = gp.Model()
        
        x = model.addVars([a for a in range(nb_arcs)], \
                    lb = 0.0,\
                    ub = 1.0,\
                    obj = 0.0,\
                    vtype = GRB.BINARY,\
                    name = "x") 
        
        y = model.addVars([a for a in range(nb_arcs)], \
                lb = 0.0,\
                ub = GRB.INFINITY,\
                obj = dist_arcs,\
                vtype = GRB.CONTINUOUS,\
                name = "y") 
        
        z = model.addVars([n for n in range(self.nb_nodes)], \
            lb = 0.0,\
            ub = 1.0,\
            obj = [median_house_cost for n in range(self.nb_nodes)],\
            vtype = GRB.BINARY,\
            name = "z") 
        
        w = model.addVars([n for n in range(self.nb_nodes)], \
            lb = 0.0,\
            ub = GRB.INFINITY,\
            obj = 0.0,\
            vtype = GRB.CONTINUOUS,\
            name = "w") 
        
        model.addConstrs(gp.quicksum(y[a] for a in range(nb_arcs) if arcs[a][0] == n ) \
                         - gp.quicksum(y[a] for a in range(nb_arcs) if arcs[a][1] == n ) \
                         == demand[n] - w[n] for n in range(self.nb_nodes))
        
        model.addConstrs(w[n] <= max_node_capacity*z[n] for n in range(self.nb_nodes))
        
        model.addConstrs(y[a] <= max_arc_capacity*x[a] for a in range(nb_arcs))
        
        model.addConstr(z.sum("*") == self.nb_facilities)
        
        model.optimize()
        
        self.isFacility = list()
        if model.Status == GRB.OPTIMAL:
            for n in range(self.nb_nodes):
                self.isFacility.append(z[n].x)
        else:
            print("Problem, the model has not obtained the optimal solution.")
    
    def writeFile(self):
        # Writing file
        name = '../instances/usa/' + str(self.nb_nodes) + 'nodes' + str(self.nb_facilities) + 'facilities.txt'
        with open(name, 'w') as f:
            f.write( "%d %d %d %d %d\n" % (self.nb_nodes,self.nb_edges,self.nb_facilities,nb_disr_type,nb_disr_level))
            for i in range(self.nb_nodes):
                f.write("%s %d %d %d %d %d\n" % (self.cities['city'].iloc[i].replace(" ", ""),round(self.cities['population'].iloc[i]/10000),self.disr1[i],self.disr2[i],self.disr3[i],self.isFacility[i])) 
            for a,b in self.edges:
                dist =  math.pow((self.points[a,0] - self.points[b,0])**2 + (self.points[a,1] - self.points[b,1])**2, 1/2)/1000
                f.write("%s %s %.3f\n" % (self.cities['city'].iloc[a].replace(" ", ""),self.cities['city'].iloc[b].replace(" ", ""),dist))
            for d in range(nb_disr_type):
                f.write("%s %.3f\n" % (self.cities_and_disr_places['city'].iloc[self.nb_nodes+d].replace(" ", ""),prob_disr[d]))

def main():
    arg = parse_arg()
    data = Data(arg)
    
if __name__ == "__main__":
    main()

