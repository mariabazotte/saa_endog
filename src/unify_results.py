import pandas as pd
import openpyxl
import argparse
import numpy as np
import math

########## Paramters configuration #########
seeds         = [0,1000,2000,3000,4000]                                                 # Seeds
network_nodes = ["15nodes","25nodes","30nodes","39nodes","48nodes"]                     # Nodes
facilities    = ["4facilities","5facilities"]                                           # Facilities
networks      = [[node,facility] for node in network_nodes for facility in facilities]  # Networks

parameters    = ["params"]                    # Parameters file
budgets       = [0.5]                         # Budgets
capas         = [0,2,3,4]                     # Capacity states (0 is for the normal)
val_distr     = [0,1,2,3]                     # Distributions
nb_problems         = 50                      # Number of problems
chosen_nb_scenarios = 750                     # Chosen number of scenarios (between list of scenarios) for a detailed analysis
nb_val_scenarios    = 50000                   # Number of scenarios for validation
edgemult = 10                                 # Cost multiplier for edge 

########## Folders #########
folder_distr  = ["discreteselection/","binomial/","normal/","stdnormalization/"]  # Folders for distributions
linearization_folder = "discreteselection/"

############ Instances and Distributions information ############
distributions = ["Discrete selection","Binomial","Normal","Std. Normalization"]
nb_nodes      = {"10nodes":10,"15nodes":15,"25nodes":25,"30nodes":30,"39nodes":39,"48nodes":48}
nb_edges      = {"10nodes":22,"15nodes":36,"25nodes":66,"30nodes":80,"39nodes":106,"48nodes":131}
nb_facilities = {"4facilities":4,"5facilities":5}

############ Name of file #############
excel_file    = 'results.xlsx'
      
def getDistributions(capa):
    if(capa == 0):
        if(2 in val_distr):
            return [2]
        else:
            return []
    else:
        l = []
        if(0 in val_distr): 
            l.append(0)
        if(1 in val_distr):
            l.append(1)
        if(3 in val_distr):
            l.append(3)
        return l

dict_tables = {"seed":[],"network_file":[],"parameters_file":[],"nb_nodes":[],"nb_edges":[],"nb_facilities":[],\
               "nb_capa_states":[],"budget":[],"nb_scenarios":[],"lb":[],"ub":[],"gap":[],"nbnodes":[],"time":[],\
               "saalb":[],"saaub":[],"abssaagap":[],"saagap":[],"stdsaalb":[],"stdsaaub":[],"stdsaagap":[],"saastagap":[],\
               "saanbnodes":[],"stdsaanbnodes":[],"saatime":[],"gaplb":[],"gapub":[]}
dict_mean_tables = {"network_file":[],"parameters_file":[],"nb_nodes":[],"nb_edges":[],"nb_facilities":[],"nb_capa_states":[],\
                    "budget":[],"nb_scenarios":[],"lb":[],"intlb":[],"ub":[],"intub":[],"gap":[],"intgap":[],"nbnodes":[],"intnbnodes":[],\
                    "time":[],"inttime":[],"saalb":[],"intsaalb":[],"saaub":[],"intsaaub":[],"abssaagap":[],"intabssaagap":[],\
                    "saagap":[],"intsaagap":[],"stdsaalb":[],"intstdsaalb":[],"stdsaaub":[],"intstdsaaub":[],"stdsaagap":[],"intstdsaagap":[],\
                    "saastagap":[],"intsaastagap":[],"saanbnodes":[],"intsaanbnodes":[],"stdsaanbnodes":[],"intstdsaanbnodes":[],"saatime":[],"intsaatime":[],\
                    "gaplb":[],"intgaplb":[],"gapub":[],"intgapub":[]}

dict_tables_ev = {"seed":[],"distribution":[],"network_file":[],"parameters_file":[],"nb_nodes":[],"nb_edges":[],"nb_facilities":[],"nb_capa_states":[],\
                  "budget":[],"nb_scenarios":[],"saalb":[],"saaub":[],"abssaagap":[],"saagap":[],"stdsaalb":[],"stdsaaub":[],"stdsaagap":[],"saastagap":[],\
                  "saanbnodes":[],"stdsaanbnodes":[],"nbsaanotsolved":[],"saatime":[],"EV":[],"EEV":[],"VSS":[],"EEV_SAA":[],"STD_EEV_SAA":[],"VSS_SAA":[],"STD_VSS_SAA":[]}
dict_mean_tables_ev = {"distribution":[],"network_file":[],"parameters_file":[],"nb_nodes":[],"nb_edges":[],"nb_facilities":[],"nb_capa_states":[],\
                       "budget":[],"nb_scenarios":[],"saalb":[],"intsaalb":[],"saaub":[],"intsaaub":[],"abssaagap":[],"intabssaagap":[],"saagap":[],"intsaagap":[],\
                       "stdsaalb":[],"intstdsaalb":[],"stdsaaub":[],"intstdsaaub":[],"stdsaagap":[],"intstdsaagap":[],"saastagap":[],"intsaastagap":[],\
                       "saanbnodes":[],"intsaanbnodes":[],"stdsaanbnodes":[],"intstdsaanbnodes":[],"nbsaanotsolved":[],"intnbsaanotsolved":[],"saatime":[],"intsaatime":[],"EV":[],"intEV":[],\
                       "EEV":[],"intEEV":[],"VSS":[],"intVSS":[],"EEV_SAA":[],"intEEV_SAA":[],"STD_EEV_SAA":[],"intSTD_EEV_SAA":[],\
                       "VSS_SAA":[],"intVSS_SAA":[],"STD_VSS_SAA":[],"intSTD_VSS_SAA":[]}

dict_tables_comp_vi = {"seed":[],"distribution":[],"network_file":[],"parameters_file":[],"nb_nodes":[],"nb_edges":[],"nb_facilities":[],"nb_capa_states":[],\
                       "budget":[],"1branchnode":[],"1stdbranchnode":[],"1saagap":[],"1stdsaagap":[],"1saastagap":[],"1nbsaanotsolved":[],"1time":[],\
                       "0branchnode":[],"0stdbranchnode":[],"0saagap":[],"0stdsaagap":[],"0saastagap":[],"0nbsaanotsolved":[],"0time":[]}
dict_mean_tables_comp_vi = {"distribution":[],"network_file":[],"parameters_file":[],"nb_nodes":[],"nb_edges":[],"nb_facilities":[],"nb_capa_states":[],\
                            "budget":[],"1branchnode":[],"1intbranchnode":[],"1stdbranchnode":[],"1intstdbranchnode":[],"1saagap":[],"1intsaagap":[],\
                            "1stdsaagap":[],"1intstdsaagap":[],"1saastagap":[],"1intsaastagap":[],"1nbsaanotsolved":[],"1intnbsaanotsolved":[],"1time":[],"1inttime":[],\
                            "0branchnode":[],"0intbranchnode":[],"0stdbranchnode":[],"0intstdbranchnode":[],"0saagap":[],"0intsaagap":[],\
                            "0stdsaagap":[],"0intstdsaagap":[],"0saastagap":[],"0intsaastagap":[],"0nbsaanotsolved":[],"0intnbsaanotsolved":[],"0time":[],"0inttime":[]}

# Considering 95% confidence
def gettstudent():
    if(len(seeds)==1):
        return 0
    if(len(seeds)==2):
        return 12.71
    if(len(seeds)==3):
        return 4.30
    if(len(seeds)==4):
        return 3.18
    if(len(seeds)==5):
        return 2.78
    if(len(seeds)==6):
        return 2.57
    if(len(seeds)==7):
        return 2.45
    if(len(seeds)==8):
        return 2.36
    if(len(seeds)==7):
        return 2.31
    if(len(seeds)==10):
        return 2.26
    return 0

def setdict(seed,network_pair,parameters_file,capa,budget):
    dict_tables["seed"].append(seed)
    dict_tables["network_file"].append(network_pair[0]+network_pair[1])
    dict_tables["parameters_file"].append(parameters_file)
    dict_tables["nb_nodes"].append(nb_nodes[network_pair[0]])
    dict_tables["nb_edges"].append(nb_edges[network_pair[0]])
    dict_tables["nb_facilities"].append(network_pair[1])
    dict_tables["nb_capa_states"].append(capa)
    dict_tables["budget"].append(budget)

def setmeandict(network_pair,parameters_file,capa,budget):
    dict_mean_tables["network_file"].append(network_pair[0]+network_pair[1])
    dict_mean_tables["parameters_file"].append(parameters_file)
    dict_mean_tables["nb_nodes"].append(nb_nodes[network_pair[0]])
    dict_mean_tables["nb_edges"].append(nb_edges[network_pair[0]])
    dict_mean_tables["nb_facilities"].append(network_pair[1])
    dict_mean_tables["nb_capa_states"].append(capa)
    dict_mean_tables["budget"].append(budget)
    
    dict_mean_tables["nb_scenarios"].append(dict_tables["nb_scenarios"][-1])
    for opt in ["lb","ub","gap","nbnodes","time","saalb","saaub","abssaagap","saagap","stdsaalb","stdsaaub","stdsaagap","saastagap","saanbnodes","stdsaanbnodes","saatime","gaplb","gapub"]:
        dict_mean_tables[opt].append(np.mean(dict_tables[opt][-len(seeds):]))
        dict_mean_tables["int" + opt].append( str(round(dict_mean_tables[opt][-1],1)) + "+/-" + str(round((gettstudent()*np.std(dict_tables[opt][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),1)))
    
def setlinearization(nb_sce,lb,ub,gap,nbnodes,time):
    dict_tables["nb_scenarios"].append(nb_sce)
    dict_tables["lb"].append(lb)
    dict_tables["ub"].append(ub)
    dict_tables["gap"].append(gap)
    dict_tables["nbnodes"].append(nbnodes)
    dict_tables["time"].append(time)

def setsaa(saalb,saaub,saagap,stdsaalb,stdsaaub,stdsaagap,statgap,saanbnodes,stdsaanbnodes,saatime):
    dict_tables["saalb"].append(saalb)
    dict_tables["saaub"].append(saaub)
    dict_tables["saagap"].append(saagap)
    dict_tables["abssaagap"].append(saagap*saaub/100)
    dict_tables["stdsaalb"].append(stdsaalb)
    dict_tables["stdsaaub"].append(stdsaaub)
    dict_tables["stdsaagap"].append(stdsaagap)
    dict_tables["saastagap"].append(statgap)
    dict_tables["saanbnodes"].append(saanbnodes)
    dict_tables["stdsaanbnodes"].append(stdsaanbnodes)
    dict_tables["saatime"].append(saatime)
    
    size = len(dict_tables["saalb"])
    dict_tables["gaplb"].append(100*(dict_tables["lb"][size-1] - dict_tables["saalb"][size-1])/dict_tables["lb"][size-1])
    dict_tables["gapub"].append(100*(dict_tables["ub"][size-1] - dict_tables["saaub"][size-1])/dict_tables["ub"][size-1])

def setevdict(distribution,seed,network_pair,parameters_file,capa,budget):
    dict_tables_ev["distribution"].append(distribution)
    dict_tables_ev["seed"].append(seed)
    dict_tables_ev["network_file"].append(network_pair[0]+network_pair[1])
    dict_tables_ev["parameters_file"].append(parameters_file)
    dict_tables_ev["nb_nodes"].append(nb_nodes[network_pair[0]])
    dict_tables_ev["nb_edges"].append(nb_edges[network_pair[0]])
    dict_tables_ev["nb_facilities"].append(network_pair[1])
    dict_tables_ev["nb_capa_states"].append(capa)
    dict_tables_ev["budget"].append(budget)
    
def setevmeandict(distribution,network_pair,parameters_file,capa,budget):
    dict_mean_tables_ev["distribution"].append(distribution)
    dict_mean_tables_ev["network_file"].append(network_pair[0]+network_pair[1])
    dict_mean_tables_ev["parameters_file"].append(parameters_file)
    dict_mean_tables_ev["nb_nodes"].append(nb_nodes[network_pair[0]])
    dict_mean_tables_ev["nb_edges"].append(nb_edges[network_pair[0]])
    dict_mean_tables_ev["nb_facilities"].append(network_pair[1])
    dict_mean_tables_ev["nb_capa_states"].append(capa)
    dict_mean_tables_ev["budget"].append(budget)
    
    dict_mean_tables_ev["nb_scenarios"].append(dict_tables_ev["nb_scenarios"][-1])
    for opt in ["saalb","saaub","abssaagap","saagap","stdsaalb","stdsaaub","stdsaagap","saastagap","saanbnodes","stdsaanbnodes","nbsaanotsolved","saatime","EV","EEV","VSS","EEV_SAA","STD_EEV_SAA","VSS_SAA","STD_VSS_SAA"]:
        dict_mean_tables_ev[opt].append(np.mean(dict_tables_ev[opt][-len(seeds):]))
        dict_mean_tables_ev["int" + opt].append(str(round(dict_mean_tables_ev[opt][-1],1)) +"+/-" + str(round((gettstudent()*np.std(dict_tables_ev[opt][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),1)))
    
def setev(nb_scenarios,ev,eev,vss,eevsaa,vareevsaa,vsssaa,stdvsssaa):
    dict_tables_ev["nb_scenarios"].append(nb_scenarios)
    dict_tables_ev["EV"].append(ev)
    dict_tables_ev["EEV"].append(eev)
    dict_tables_ev["VSS"].append(vss)
    dict_tables_ev["EEV_SAA"].append(eevsaa)
    dict_tables_ev["STD_EEV_SAA"].append(vareevsaa)
    dict_tables_ev["VSS_SAA"].append(vsssaa)
    dict_tables_ev["STD_VSS_SAA"].append(stdvsssaa)

def setevsaa(saalb,saaub,saagap,stdsaalb,stdsaaub,stdsaagap,saastagap,saanbnodes,stdsaanbnodes,nbsaanotsolved,saatime):
    dict_tables_ev["saalb"].append(saalb)
    dict_tables_ev["saaub"].append(saaub)
    dict_tables_ev["saagap"].append(saagap)
    dict_tables_ev["abssaagap"].append(saagap*saaub/100)
    dict_tables_ev["stdsaalb"].append(stdsaalb)
    dict_tables_ev["stdsaaub"].append(stdsaaub)
    dict_tables_ev["stdsaagap"].append(stdsaagap)
    dict_tables_ev["saastagap"].append(saastagap)
    dict_tables_ev["saanbnodes"].append(saanbnodes)
    dict_tables_ev["stdsaanbnodes"].append(stdsaanbnodes)
    dict_tables_ev["nbsaanotsolved"].append(nbsaanotsolved)
    dict_tables_ev["saatime"].append(saatime)

def setdictcompvi(distribution,seed,network_pair,parameters_file,capa,budget):
    dict_tables_comp_vi["distribution"].append(distribution)
    dict_tables_comp_vi["seed"].append(seed)
    dict_tables_comp_vi["network_file"].append(network_pair[0]+network_pair[1])
    dict_tables_comp_vi["parameters_file"].append(parameters_file)
    dict_tables_comp_vi["nb_nodes"].append(nb_nodes[network_pair[0]])
    dict_tables_comp_vi["nb_edges"].append(nb_edges[network_pair[0]])
    dict_tables_comp_vi["nb_facilities"].append(network_pair[1])
    dict_tables_comp_vi["nb_capa_states"].append(capa)
    dict_tables_comp_vi["budget"].append(budget)

def setsaacompvi(saagap,stdsaagap,branchnode,stdbranchnode,saastagap,nbsaanotsolved,saatime,vi):
    dict_tables_comp_vi[str(vi)+"saagap"].append(saagap)
    dict_tables_comp_vi[str(vi)+"stdsaagap"].append(stdsaagap)
    dict_tables_comp_vi[str(vi)+"branchnode"].append(branchnode)
    dict_tables_comp_vi[str(vi)+"stdbranchnode"].append(stdbranchnode)
    dict_tables_comp_vi[str(vi)+"saastagap"].append(saastagap)
    dict_tables_comp_vi[str(vi)+"nbsaanotsolved"].append(nbsaanotsolved)
    dict_tables_comp_vi[str(vi)+"time"].append(saatime)

def setmeancompvi(distribution,network_pair,parameters_file,capa,budget):
    dict_mean_tables_comp_vi["distribution"].append(distribution)
    dict_mean_tables_comp_vi["network_file"].append(network_pair[0]+network_pair[1])
    dict_mean_tables_comp_vi["parameters_file"].append(parameters_file)
    dict_mean_tables_comp_vi["nb_nodes"].append(nb_nodes[network_pair[0]])
    dict_mean_tables_comp_vi["nb_edges"].append(nb_edges[network_pair[0]])
    dict_mean_tables_comp_vi["nb_facilities"].append(network_pair[1])
    dict_mean_tables_comp_vi["nb_capa_states"].append(capa)
    dict_mean_tables_comp_vi["budget"].append(budget)
    
    for opt in ["saagap","stdsaagap","nbsaanotsolved","time","branchnode","stdbranchnode","saastagap"]:
        for vi in [0,1]:
            dict_mean_tables_comp_vi[str(vi) + opt].append(np.mean(dict_tables_comp_vi[str(vi) + opt][-len(seeds):]))
            dict_mean_tables_comp_vi[str(vi) + "int" + opt].append( str(round(dict_mean_tables_comp_vi[str(vi) + opt][-1],1)) + "+/-" + str(round((gettstudent()*np.std(dict_tables_comp_vi[str(vi) + opt][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),1)))
    
def main():
    path = '../results/compiled/' + str(nb_problems) + '_' + str(chosen_nb_scenarios) + "_" + str(nb_val_scenarios) + "_" + excel_file
    with pd.ExcelWriter(path) as writer:
        for network_pair in networks:
            network = network_pair[0] + network_pair[1]
            network_df = pd.DataFrame() 
            for parameter in parameters:
                for budget in budgets:
                    for capa in capas:
                        if (capa != 0):
                            for seed in seeds:
                                setdict(seed,network_pair,parameter,capa,budget)
                                linearization_file = "../results/linearization/" + linearization_folder + network + "_" + parameter + "_budget" + str(int(budget*100)) + "_capa" + str(capa) + "_edgemult" + str(edgemult)
                    
                                name_file = linearization_file + "_call0_ppd1_ev1_seed" + str(seed) + ".csv" # L-Shaped linearization (without callback)
                                df = pd.read_csv(name_file,delimiter=';')
                                
                                name_file = linearization_file + "_call1_ppd1_ev1_seed" + str(seed) + ".csv" # L-Shaped linearization (with callback)
                                # df = pd.read_csv(name_file,delimiter=';')
                                
                                network_df = pd.concat([network_df,df]) 
                                setlinearization(df["TOTAL_NB_SCE"][0],\
                                                float(df["LB"][0]),\
                                                float(df["UB"][0]),\
                                                float(df["GAP"][0])*100,\
                                                float(df["NB_BRANCH_NODES"][0]),\
                                                float(df["TIME"][0]))

                        distr = getDistributions(capa) 
                        for d in distr: 
                            expected_file = "../results/" + "expected/" + folder_distr[d] + network + "_" + parameter + "_budget" + str(int(budget*100)) + "_capa" + str(capa) + "_edgemult" + str(edgemult) +"_"
                            saa_file      = "../results/" + "saa/"      + folder_distr[d] + network + "_" + parameter + "_budget" + str(int(budget*100)) + "_capa" + str(capa) + "_edgemult" + str(edgemult) +"_"
                            for seed in seeds:
                                setevdict(distributions[d],seed,network_pair,parameter,capa,budget)
                                name_file = expected_file + "seed" + str(seed) + ".csv" # Expected 
                                mean_df = pd.read_csv(name_file,delimiter=';')
                                network_df = pd.concat([network_df,mean_df])
                            
                                total_nb_scenarios = pow(capa+1,nb_facilities[network_pair[1]])*4
                            
                                if(d == 1 or d == 3):
                                    setdictcompvi(distributions[d],seed,network_pair,parameter,capa,budget)
                                    for vi in [0,1]:
                                        name_file = saa_file + "vi" + str(vi) + "_"  + str(nb_problems) + "_" + str(chosen_nb_scenarios) + "_" + str(nb_val_scenarios) + "_seed" + str(seed) + ".csv" # SAA
                                        df = pd.read_csv(name_file,delimiter=';')
                                        network_df = pd.concat([network_df,df]) 
                                        
                                        if(vi == 0):
                                            setevsaa(float(df["LB"][0]),\
                                                    float(df["UB"][0]),\
                                                    float(df["GAP"][0])*100,\
                                                    math.sqrt(float(df["VAR_LB"][0])),\
                                                    math.sqrt(float(df["VAR_UB"][0])),\
                                                    math.sqrt(float(df["VAR_GAP"][0])),\
                                                    float(df["STAT_GAP"][0]),\
                                                    float(df["NB_BRANCH_NODES"][0]),\
                                                    math.sqrt(float(df["VAR_BRANCH_NODES"][0])),\
                                                    float(df["NB_SAA_NOT_SOLVED"][0]),\
                                                    float(df["TIME"][0]))  

                                        setsaacompvi(float(df["GAP"][0])*100,\
                                                    math.sqrt(float(df["VAR_GAP"][0])),\
                                                    float(df["NB_BRANCH_NODES"][0]),\
                                                    math.sqrt(float(df["VAR_BRANCH_NODES"][0])),\
                                                    float(df["STAT_GAP"][0]),\
                                                    float(df["NB_SAA_NOT_SOLVED"][0]),\
                                                    float(df["TIME"][0]),vi)                     
                                else:
                                    name_file = saa_file + str(nb_problems) + "_" + str(chosen_nb_scenarios) + "_" + str(nb_val_scenarios) + "_seed" + str(seed) + ".csv" # SAA
                                    df = pd.read_csv(name_file,delimiter=';')
                                    network_df = pd.concat([network_df,df])
                                    
                                    if(d == 0):
                                        setsaa(float(df["LB"][0]),\
                                            float(df["UB"][0]),\
                                            float(df["GAP"][0])*100,\
                                            math.sqrt(float(df["VAR_LB"][0])),\
                                            math.sqrt(float(df["VAR_UB"][0])),\
                                            math.sqrt(float(df["VAR_GAP"][0])),\
                                            float(df["STAT_GAP"][0]),\
                                            float(df["NB_BRANCH_NODES"][0]),\
                                            math.sqrt(float(df["VAR_BRANCH_NODES"][0])),\
                                            float(df["TIME"][0]))
                                        
                                    setevsaa(df["LB"][0],\
                                            df["UB"][0],\
                                            float(df["GAP"][0])*100,\
                                            math.sqrt(df["VAR_LB"][0]),\
                                            math.sqrt(df["VAR_UB"][0]),\
                                            math.sqrt(df["VAR_GAP"][0]),\
                                            float(df["STAT_GAP"][0]),\
                                            float(df["NB_BRANCH_NODES"][0]),\
                                            math.sqrt(float(df["VAR_BRANCH_NODES"][0])),\
                                            float(df["NB_SAA_NOT_SOLVED"][0]),\
                                            float(df["TIME"][0]))
                                    
                                if((d == 0 or d == 1 or d == 3) and total_nb_scenarios < 20000):
                                    setev(mean_df["TOTAL_NB_SCE"][2],\
                                          mean_df["UB"][0],\
                                          mean_df["UB"][2],\
                                          100*(mean_df["UB"][2]-dict_tables_ev["saaub"][-1])/mean_df["UB"][2],\
                                          mean_df["UB"][1],\
                                          math.sqrt(mean_df["VAR_UB"][1]),\
                                          100*(mean_df["UB"][1]-dict_tables_ev["saaub"][-1])/mean_df["UB"][1],\
                                          math.sqrt(mean_df["VAR_UB"][1]+0))
                                else:
                                    setev(0,mean_df["UB"][0],0,0,\
                                          mean_df["UB"][1],\
                                          math.sqrt(mean_df["VAR_UB"][1]),\
                                          100*(mean_df["UB"][1]-dict_tables_ev["saaub"][-1])/mean_df["UB"][1],\
                                          math.sqrt(mean_df["VAR_UB"][1]+(dict_tables_ev["stdsaaub"][-1])**2))
                                    
                            setevmeandict(distributions[d],network_pair,parameter,capa,budget)
                            if(d == 1 or d == 3):
                                setmeancompvi(distributions[d],network_pair,parameter,capa,budget)
                        if(capa != 0):
                            setmeandict(network_pair,parameter,capa,budget) 
            sheet = network
            network_df.to_excel(writer, sheet_name=sheet) 
            
        # For comparing discrete choice saa with baseline
        tables_df = pd.DataFrame(dict_tables)  
        sheet = "tables_per_seed"
        tables_df.to_excel(writer,sheet_name=sheet)
        if(len(seeds)>1):
            mean_table_df = pd.DataFrame(dict_mean_tables)
            sheet = "tables_mean"
            mean_table_df.to_excel(writer,sheet_name=sheet)
                
        # For comparing with EV and EEV
        ev_tables_df = pd.DataFrame(dict_tables_ev) 
        sheet = "ev_tables_per_seed"
        ev_tables_df.to_excel(writer,sheet_name=sheet)
        if(len(seeds)>1):
            mean_ev_tables_df = pd.DataFrame(dict_mean_tables_ev) 
            sheet = "ev_tables_mean"
            mean_ev_tables_df.to_excel(writer,sheet_name=sheet)
            
        # For analysing the valid inequality
        vi_tables_df = pd.DataFrame(dict_tables_comp_vi)
        sheet = "vi_tables_per_seed"
        vi_tables_df.to_excel(writer,sheet_name=sheet)
        if(len(seeds)>1):
            mean_vi_tables_df = pd.DataFrame(dict_mean_tables_comp_vi) 
            sheet = "vi_tables_mean"
            mean_vi_tables_df.to_excel(writer,sheet_name=sheet)
    
if __name__ == "__main__":
    main()
