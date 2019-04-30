import os
import numpy as np
import numpy.testing as npt
import glob
import re
import argparse
import vtk
import sys
import shutil
import flopy.utils.binaryfile as fpu
import pandas as pd
import flopy


from math import sin, pi
    
# cd : critical duration  chosen to be 7days
# h : d
# hc : crtical value for the head of water table
# ratio = 1 : deltaDc / cd = ratio 
def costDuration(d, dc, ratio):
    deltaDc = dc * ratio
    if (d >= (dc-(deltaDc/2))) and (d <= (dc+(deltaDc/2))):
        expression = pi * ((d- (dc - (deltaDc/2)) )/deltaDc)
        return sin(expression)
    else :
        return 0.0 



#repo = "/DATA/These/Projects/Model/docker/app/Agon-Coutainville/model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodtenyear_Step1/"
#modelname = "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodtenyear_Step1"
upperLimitForFloodZone = 0.3




# os.chdir(repo)
# vtuFiles = glob.glob("output_files/VTU_Grid.vtu")

# if (len(vtuFiles) == 0):
#     print("There is no VTU file in the repository of the given modelname. Check if the files have indeed been generated in the previous phase of the pipeline.")
#     sys.exit(0)
    
#     # One file stores a watertable of the geographical zone (matrix) for a stress period
    

#     #shutil.copyfile(vtuFiles[0], "VTU_WaterTable_FloodDuration.vtu")
# #fin = open(repo + "output_files/VTU_Grid.vtu", "r" )
# fin = open (repo + vtuFiles[0], "r")
# data_list = fin.readlines()
# fin.close()

# index = [x for x in range(len(data_list)) if '</celldata>' in data_list[x].lower()]

# del data_list[4:index[0]]

# fout = open("VTU_Grid_" + modelname + "_FloodDuration_" + str(upperLimitForFloodZone) + ".vtu", "w")
# fout.writelines(data_list)
# fout.close()

    
def getFloodDuration(repo, modelname):  
    os.chdir(repo)
    hds = fpu.HeadFile(modelname + '.hds')
    times = hds.get_times()
    kstpkper = hds.get_kstpkper()


    mf = flopy.modflow.Modflow.load(repo + modelname + '.nam')
    dis = flopy.modflow.ModflowDis.load(repo + modelname + '.dis', mf)
    topo = dis.top._array


    sites = pd.read_csv("/DATA/These/Projects/Model/docker/app/" + "data/" + "study_sites.txt", sep='\s+', header=0, index_col=0)
    site_number = 2 #Select site number
    port=int(sites._get_values[site_number,5])
    ram = pd.read_csv("/DATA/These/Projects/Model/docker/app/"+"data/RAM.csv", delimiter=";", header=0)
    sea_level = int(ram.NM_IGN[port-1]) # 0 ici


    floods = []

    for period in range(1, len(kstpkper)): # Pas initialisation
        
        
        # si cellule sèche, value couche 2
        
        head = hds.get_data(kstpkper=kstpkper[period])
        couches = head.shape[0]
        lines = head.shape[1]
        cols = head.shape[2]
        flood = np.zeros((lines, cols))
        for l in range(lines):
            for c in range(cols):
                ztopo = topo[l][c]
                h = head[0][l][c]
                try:
                    npt.assert_approx_equal(h,-1e+30)
                except:
                    try:
                        npt.assert_approx_equal(h,sea_level)
                        print("sea zone")
                    except:
                        if ((ztopo - h) <= upperLimitForFloodZone):
                            #print("flood")
                            if (period == 0):
                                flood[l][c] = round(times[period])
                            else :
                                flood[l][c] = round(times[period] - times[period-1])                        
                    
                        #else:
                            #print("no flood")
        #floods[period] = np.array()
        floods.append(flood)
        #print(floods[period-1])
        #print(flood)
    floods = np.array(floods)
    floodstot = np.sum(floods, axis=0)
    return floodstot



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-repo", "--repo", type=str)
    parser.add_argument("-m", "--modelname", type=str, required=True)
    parser.add_argument("-repoSimu", "--repoSimu", type=str)
    parser.add_argument("-ref", "--refsimu", type=str, required=True)
   
    args = parser.parse_args()

    repo = args.repo
    modelname = args.modelname
    refSimu = args.refsimu
    repoSimu = args.repoSimu

    floodRef = getFloodDuration(repoSimu, refSimu)
    floodAlt = getFloodDuration(repo, modelname)

    mae = 0
    mre = 0
    rmse = 0

    for l in range(floodRef.shape[0]):
        for c in range(floodRef.shape[1]):
            diff = floodRef[l][c] - floodAlt[l][c]
        mae += abs(diff)
        mre += abs(diff / max(1, (floodRef[l][c] - floodAlt[l][c])/2))
        rmse += diff**2

    mae = mae / tailleRef
    rmse = math.sqrt(rmse / tailleRef)

    print ("mae" + str(mae))
    print ("mre" + str(mre))
    print ("rmse" + str(rmse))
   