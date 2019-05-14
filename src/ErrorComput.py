
"""
    Provides functions to compute error rate between the reference simulation and an alternative simulation.
"""

import os
import numpy as np
import re
import argparse
import sys
import flopy.utils.binaryfile as fpu
import math
import csv
import pandas as pd
import math
import flopy
import pickle
import statistics

__author__ = "June Sallou"
__maintainer__ = "June Sallou"
__credits__ = ["June Sallou"]
__license__ = "MIT"
__version__ = "0.0.3"
__date__ = "04/30/2019"
__email__ = "june.benvegnu-sallou@univ-rennes1.fr"


mainAppRepo = os.getcwd() + "/"


def getPathToSimulationDirectoryFromModelname(modelname, site_name):
    return os.path.join(mainAppRepo, site_name, modelname)

def getSiteNameFromSiteNumber(site_number):
    sites = pd.read_csv(mainAppRepo + 'data/study_sites.txt', sep='\\s+', header=0, index_col=0)
    site_name = sites.index._data[site_number] + '/'
    return site_name

        
def getNonDryCellHdsValue(hds, nrow, ncol, nlayer):
    layer = 0
    h = hds[layer][nrow][ncol]
    while (math.isclose(abs(h)/1e+30, 1, rel_tol=1e-3)) and layer<nlayer:
        if layer == nlayer-1:
            print("cell completely dry")
        else:
            h = hds[layer+1][nrow][ncol]
            layer+=1
    return h



def getWeightToSurface(zs, h, dc, alpha):
    """
        zs : value of soil surface at the point x
        h : head value for watertable at the x point
        dc : critical depth
        alpha : ratio of critical depth for calculating the width of transition zone
    """
    ddc = alpha * dc #Alpha must be not null
    
    borneInf = zs - (dc + (ddc / 2))
    borneSup =  zs - (dc - (ddc / 2))
    
    if h <= borneInf:
        Ws = 0
    elif h >= borneSup:
        Ws = 1
    else :
        Ws = math.sin((math.pi * (h - borneSup)) / (2*ddc)) 
    
    return Ws


def getSoilSurfaceValuesForASimulation(repoSimu, modelname):
    """
        Retrieve the matrix of the topographical altitude.
    """
    mf = flopy.modflow.Modflow.load(repoSimu + "/" + modelname + '.nam')
    dis = flopy.modflow.ModflowDis.load(repoSimu + "/" + modelname + '.dis', mf)
    topo = dis.top._array

    return topo


def computeFloodsAndHErrorFromModelnamesByInterpolation(ref, modelname, site_number,timestep=1):
    
    site_name = getSiteNameFromSiteNumber(site_number)
    
    # Get data for reference simulation
    ## Path to repository
    repoRef = getPathToSimulationDirectoryFromModelname(ref, site_name) 
    ## Topographical altitude
    topoRef = getSoilSurfaceValuesForASimulation(repoRef, ref)
    ## Watertable altitude
    refHds = fpu.HeadFile(repoRef + '/' + ref + '.hds')

    # Get data for alternate simulation
    ## Path to repository
    repoSimu = getPathToSimulationDirectoryFromModelname(modelname, site_name)
    ## Topographical altitude
    topoSimu = getSoilSurfaceValuesForASimulation(repoSimu, modelname)
    # Get heads values for simulation
    simuHds = fpu.HeadFile(repoSimu + '/' + modelname + '.hds')
    ## Duration of simulated periods 
    simuTimes = simuHds.get_times()
    
    # Comput parameters
    startTime = 0
    endTime = 15340
    dc = 0.3
    alpha = float(1/3)
    
    # Initialisation of error computing variables
    floods = {}
    smWs = 0
    sherrorsup = 0


    #For every day, we have to compare the reference and alternate simulations
    for day in range(startTime, endTime+1):
        # print(day)
    
        # Retrieve the watertable altitude values for the reference simulation for the considered day 
        refHead = refHds.get_data(kstpkper=(0, day))

        # Compute the number of the simulated period of the alternate simulation for the considered day
        nbPeriod = 0
        while (simuTimes[nbPeriod] < day+1) and (nbPeriod < len(simuTimes)):
            nbPeriod+=1
        # print("nbPeriod : " + str(nbPeriod))
        
        # When the considered day match with the duration of the corresponding simulated period
        # We retrieve the value of the watertable altitude from the corresponding matrix
        if math.isclose(simuTimes[nbPeriod], day+1, rel_tol=1e-3): 
            altHeadSup = simuHds.get_data(kstpkper=(timestep-1, nbPeriod))       
            altHeadInf = altHeadSup
            duree = int(simuTimes[nbPeriod])
            pas = 0
        # Otherwise, we have to interpolate the watertable altitude value for the considered day
        else :
            # The considered day is situated between the simulated period number 'nbPeriod-1' and number 'nbPeriod'
            altHeadSup = simuHds.get_data(kstpkper=(timestep-1, nbPeriod))        
            altHeadInf = simuHds.get_data(kstpkper=(timestep-1, nbPeriod-1))
            duree = int(simuTimes[nbPeriod] - simuTimes[nbPeriod-1])
            pas = day - simuTimes[nbPeriod-1]

        nbrowtot = altHeadInf.shape[1]        
        nbcoltot = altHeadInf.shape[2]
        
        # We want to store the presence of cells part of a flood episode
        flood = {}
        # We go through all the cells of the matrix representing the study site
        for nrow in range(nbrowtot):
            flood[nrow] = {}
            for ncol in range(nbcoltot):
                ss = getNonDryCellHdsValue(altHeadInf, nrow, ncol, altHeadInf.shape[0]) # Watertable altitude value for simulated period with duration lower than considered day
                se = getNonDryCellHdsValue(altHeadSup, nrow, ncol, altHeadInf.shape[0]) # Watertable altitude value for simulated period with duration higher than considered day
                ajoutSimu = (se - ss) / duree
                s = ss + (ajoutSimu * pas) # Watertable altitude value for considered day being interpolated

                d = topoSimu[nrow][ncol] - s # deapth : topographical altitude  - watertable altitude
                
                # The cell is considered as undergoing a flood episode
                if d <= dc:
                    flood[nrow][ncol] = 1                    

                r= getNonDryCellHdsValue(refHead, nrow, ncol,refHead.shape[0]) # Watertable altitude value for simulated period for reference simulation
                
                WsRef = getWeightToSurface(topoRef[nrow][ncol], r, dc, alpha)
                WsSimu = getWeightToSurface(topoSimu[nrow][ncol], s, dc, alpha)

                mWs = max(WsRef, WsSimu)
                sherrorsup += (mWs * (r-s)**2)
                smWs += mWs 
        floods[day] = flood
    f = open(repoSimu + "/" + modelname + '_floods_pickle_dict.txt','wb')
    pickle.dump(floods, f)

    hErrorGlobal = math.sqrt(sherrorsup / smWs)
    
    with open(repoSimu + "/" + modelname + '_Ref_' + ref + '_errorsresult_H_light.csv', 'w') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(['H Error'])
        writer.writerow([hErrorGlobal])
    
    return nbrowtot, nbcoltot


def computeOnlyFloodsFromModelnamesByInterpolation(ref, modelname, site_number,timestep=1):
    
    site_name = getSiteNameFromSiteNumber(site_number)
    
    # Get data for reference simulation
    ## Path to repository
    # repoRef = getPathToSimulationDirectoryFromModelname(ref, site_name) 
    # ## Topographical altitude
    # topoRef = getSoilSurfaceValuesForASimulation(repoRef, ref)
    # ## Watertable altitude
    # refHds = fpu.HeadFile(repoRef + '/' + ref + '.hds')

    # Get data for alternate simulation
    ## Path to repository
    repoSimu = getPathToSimulationDirectoryFromModelname(modelname, site_name)
    ## Topographical altitude
    topoSimu = getSoilSurfaceValuesForASimulation(repoSimu, modelname)
    # Get heads values for simulation
    simuHds = fpu.HeadFile(repoSimu + '/' + modelname + '.hds')
    ## Duration of simulated periods 
    simuTimes = simuHds.get_times()
    
    # Comput parameters
    startTime = 0
    endTime = 15340
    dc = 0.3
    alpha = float(1/3)
    
    # Initialisation of error computing variables
    floods = {}
    # smWs = 0
    # sherrorsup = 0


    #For every day, we have to compare the reference and alternate simulations
    for day in range(startTime, endTime+1):
       # print(day)
    
        # Retrieve the watertable altitude values for the reference simulation for the considered day 
      #  refHead = refHds.get_data(kstpkper=(0, day))

        # Compute the number of the simulated period of the alternate simulation for the considered day
        nbPeriod = 0
        while (simuTimes[nbPeriod] < day+1) and (nbPeriod < len(simuTimes)):
            nbPeriod+=1
        # print("nbPeriod : " + str(nbPeriod))
        
        # When the considered day match with the duration of the corresponding simulated period
        # We retrieve the value of the watertable altitude from the corresponding matrix
        if math.isclose(simuTimes[nbPeriod], day+1, rel_tol=1e-3): 
            altHeadSup = simuHds.get_data(kstpkper=(timestep-1, nbPeriod))       
            altHeadInf = altHeadSup
            duree = int(simuTimes[nbPeriod])
            pas = 0
        # Otherwise, we have to interpolate the watertable altitude value for the considered day
        else :
            # The considered day is situated between the simulated period number 'nbPeriod-1' and number 'nbPeriod'
            altHeadSup = simuHds.get_data(kstpkper=(timestep-1, nbPeriod))        
            altHeadInf = simuHds.get_data(kstpkper=(timestep-1, nbPeriod-1))
            duree = int(simuTimes[nbPeriod] - simuTimes[nbPeriod-1])
            pas = day - simuTimes[nbPeriod-1]

        nbrowtot = altHeadInf.shape[1]        
        nbcoltot = altHeadInf.shape[2]
        
        # We want to store the presence of cells part of a flood episode
        flood = {}
        # We go through all the cells of the matrix representing the study site
        for nrow in range(nbrowtot):
            flood[nrow] = {}
            for ncol in range(nbcoltot):
                ss = getNonDryCellHdsValue(altHeadInf, nrow, ncol, altHeadInf.shape[0]) # Watertable altitude value for simulated period with duration lower than considered day
                se = getNonDryCellHdsValue(altHeadSup, nrow, ncol, altHeadInf.shape[0]) # Watertable altitude value for simulated period with duration higher than considered day
                ajoutSimu = (se - ss) / duree
                s = ss + (ajoutSimu * pas) # Watertable altitude value for considered day being interpolated

                d = topoSimu[nrow][ncol] - s # deapth : topographical altitude  - watertable altitude
                
                # The cell is considered as undergoing a flood episode
                if d <= dc:
                    flood[nrow][ncol] = 1
                    # print("row : ", nrow, "col: ", ncol)
                    # print("topo: ", topoSimu[nrow][ncol], "s: ", s)
                    

               # r= getNonDryCellHdsValue(refHead, nrow, ncol,refHead.shape[0]) # Watertable altitude value for simulated period for reference simulation
                
                #WsRef = getWeightToSurface(topoRef[nrow][ncol], r, dc, alpha)
               # WsSimu = getWeightToSurface(topoSimu[nrow][ncol], s, dc, alpha)

                # mWs = max(WsRef, WsSimu)
                # sherrorsup += (mWs * (r-s)**2)
                # smWs += mWs 
        floods[day] = flood
    f = open(repoSimu + "/" + modelname + '_floods_pickle_dict.txt','wb')
    pickle.dump(floods, f)

    # hErrorGlobal = math.sqrt(sherrorsup / smWs)
    
    # with open(repoSimu + "/" + modelname + '_Ref_' + ref + '_errorsresult_H_light.csv', 'w') as f:
    #     writer = csv.writer(f, delimiter=';')
    #     writer.writerow(['H Error'])
    #     writer.writerow([hErrorGlobal])
    # print("h error : " + str(hErrorGlobal))
    
    return nbrowtot, nbcoltot



def getWt1AndTsAndTrForASimulation(modelname, site_number, tsmax, nbrow, nbcol, nb_years):
    site_name = getSiteNameFromSiteNumber(int(site_number))
    repoSimu = getPathToSimulationDirectoryFromModelname(modelname, site_name)
    file = open(repoSimu + '/' + modelname + "_floods_pickle_dict.txt", 'rb')
    floods = pickle.load(file)
    file.close()
    weights = {}
    ts = {}

    
    for nrow in range(nbrow):
        weights[nrow] = {}
        ts[nrow] = {}
        for ncol in range(nbcol):
            # Indicate if a flood episode is on going during the considered day
            fSimu = False

            # Flood period
            ## Starting date
            sflood = []
            ## Ending date
            eflood = []
            mean_flood_duration = 0
            

            for day in sorted(floods): 
                if (nrow in floods[day]) and (ncol in floods[day][nrow]) and (fSimu is False):
                    # Si on entre dans une periode de saturation 
                    fSimu = True
                    sflood.append(day) # On stocke le jour de debut de la periode de saturation

                elif (nrow not in floods[day] or ncol not in floods[day][nrow]) and (fSimu):
                    # Si on sort de la periode de saturation
                    fSimu = False
                    eflood.append(day) # On stocke le jour de fin de la periode de saturation

                    
                elif (day == len(floods)-1) and (fSimu):
                    # Si dernier episode incomplet
                    eflood.append(day)

            m = 0
            while m < len(eflood):
                # Pour toutes les periodes de saturation complete (il y a une date de fin)
                # On ajoute la duree de la periode
                mean_flood_duration += (eflood[m]-sflood[m])
                #print((eflood[m]-sflood[m]))
                m+=1
            # On fait la moyenne des durees
            if len(eflood) == 0:
                mean_flood_duration = 0
            else: 
                mean_flood_duration = mean_flood_duration / nb_years 
            
            
            # Storing the values for the weights and duration of saturation periods
            wt = getWeightToSaturationDuration(tsmax, mean_flood_duration)
            if wt == 1 : # Only if ts is below 4 months per year
                weights[nrow][ncol] = wt
            

            ts[nrow][ncol] = mean_flood_duration
  
    
    fweights = open(repoSimu + "/" + modelname + '_weights_pickle_per_year.txt','wb')
    pickle.dump(weights, fweights)
    fweights.close()
    fts = open(repoSimu + "/" + modelname + '_ts_pickle_per_year.txt', 'wb') #mean_flood_duration
    pickle.dump(ts, fts)
    fts.close()




def getWeightToSaturationDuration(tsmax, mean_flood_duration):
    """
        In order to keep cells having a mean duration of flood episod inferior to 4 months (121 days) per year

    """
    if abs(mean_flood_duration) > tsmax:
        Wt = 0
    if abs(mean_flood_duration) <= tsmax:
        Wt = 1
    
    return Wt

def getWeight2ToSaturationDuration(max_diff_threshold, diff_mean_flood_duration):
    """
        In order to keep cells having a difference of mean duration btw ref and alt superior to 7 days

    """
    if abs(diff_mean_flood_duration) <= max_diff_threshold:
        Wt2 = 0
    if abs(diff_mean_flood_duration) > max_diff_threshold:
        Wt2 = 1
    
    return Wt2


def computeNumberErronousCellsForASimulationWithFiles(ref, modelname, site_number, nbrow, nbcol, tsmax, max_diff_threshold):
    site_name = getSiteNameFromSiteNumber(int(site_number))
    
    # Get data for reference simulation
    repoRef = getPathToSimulationDirectoryFromModelname(ref, site_name)
    wRef = open(repoRef + '/' + ref + "_weights_pickle_per_year.txt", 'rb') # Only weights = 1 are stored
    weightsRef = pickle.load(wRef)
    wRef.close()
    tRef = open(repoRef + '/' + ref + "_ts_pickle_per_year.txt", 'rb')
    tsRef = pickle.load(tRef)
    tRef.close()

    
    # Get data for alternative simulation
    repoSimu = getPathToSimulationDirectoryFromModelname(modelname, site_name)
    wSimu = open(repoSimu + '/' + modelname + "_weights_pickle_per_year.txt", 'rb')
    weightsSimu = pickle.load(wSimu)
    wSimu.close()
    tSimu = open(repoSimu + '/' + modelname + "_ts_pickle_per_year.txt", 'rb')
    tsSimu = pickle.load(tSimu)
    tSimu.close()

    smWt = 0
    NsErrorGlobal = 0
    NsNumerator = 0
    safe_cells_number = 0
    ref_safe_cells = 0
    simu_safe_cells = 0

    for nrow in range(nbrow):
        for ncol in range(nbcol):
            print("nb row : ",nrow, "nb col : ", ncol)
            if (ncol in weightsRef[nrow]) and (ncol in weightsSimu[nrow]): 
                mWt = 1
            elif (ncol not in weightsRef[nrow]) and (ncol not in weightsSimu[nrow]):
                mWt = 0
            else:
                mWt = 1
            print("mWt : ",mWt)

            smWt += mWt
            

            if mWt == 1:
                mean_flood_duration_diff = (tsRef[nrow][ncol] - tsSimu[nrow][ncol])
                print("mean_flood_duration_diff : ", mean_flood_duration_diff)
                wt2 = getWeight2ToSaturationDuration(max_diff_threshold, mean_flood_duration_diff)
                NsNumerator += wt2
            
            if tsRef[nrow][ncol]==0:
                print("tsRef", tsRef[nrow][ncol])
                ref_safe_cells += 1
                if tsSimu[nrow][ncol]==0:
                    safe_cells_number +=1
            if tsSimu[nrow][ncol]==0:
                simu_safe_cells +=1


    if smWt == 0:
        NsErrorGlobal = 0
        print("NsNumerator", NsNumerator, "smWt", smWt)
    else:
        print("NsNumerator", NsNumerator, "smWt", smWt)
        NsErrorGlobal = NsNumerator / smWt
        print("Number of safe cells: ", safe_cells_number)
        print("ref safe cells : ", ref_safe_cells)
        print("simu safe cells : ", simu_safe_cells)

    with open(repoSimu + "/" + modelname + '_Ref_' + ref + '_errorsresult_NS_per_year.csv', 'w') as output:
        writer = csv.writer(output, delimiter=';')
        writer.writerow(['NS Error'])
        writer.writerow([NsErrorGlobal])
    print("Ns Error value : ", NsErrorGlobal)





if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-refname", "--refname", type=str, required=False)
    parser.add_argument("-m", "--modelname", type=str, required=True)
    parser.add_argument("-site", "--sitenumber", type=int, required=True)
    parser.add_argument("-ts", "--timestep", type=int, required=False)
    parser.add_argument("-stime", "--starttime", type=int, required=False)
    parser.add_argument("-etime", "--endtime", type=int, required=False)
    parser.add_argument("-ref", "--ref", action='store_true')


    args = parser.parse_args()

    refname = args.refname
    modelname = args.modelname
    site_number = args.sitenumber
    timestep = args.timestep
    startTime = args.starttime
    endTime = args.endtime
    ref = args.ref
    
    tsmax = 121  # 4 months = 121,75 days 

    if (refname is None):
        refname = "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1"
    
    if timestep is None:
        timestep = 1

    nbrowtot = 173
    nbcoltot = 160
    nb_years = 42
    max_diff_threshold = 7


    if (ref):
        computeFloodsAndHErrorFromModelnamesByInterpolation(refname, modelname, site_number, timestep=timestep)
    else:
        getWt1AndTsAndTrForASimulation(modelname, site_number, tsmax, nbrowtot, nbcoltot, nb_years)
        computeNumberErronousCellsForASimulationWithFiles(refname, modelname, site_number, nbrowtot, nbcoltot, tsmax, max_diff_threshold)

 