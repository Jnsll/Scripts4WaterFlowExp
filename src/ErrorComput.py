#!/usr/bin/env python
"""
    Provides functions to compute error rate between the reference simulation and an alternative simulation.
"""

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
import math
import csv
import pandas as pd
import threading
import subprocess
import math
import flopy
import pickle

__author__ = "June Sallou"
__maintainer__ = "June Sallou"
__credits__ = ["June Sallou"]
__license__ = "MIT"
__version__ = "0.0.1"
__date__ = "04/30/2019"
__email__ = "june.benvegnu-sallou@univ-rennes1.fr"


mainAppRepo = os.getcwd() + "/"


def getPathToSimulationDirectoryFromModelname(modelname, site_name):
    return os.path.join(mainAppRepo, site_name, modelname)

def getSiteNameFromSiteNumber(site_number):
    sites = pd.read_csv(mainAppRepo + r"data/study_sites.txt", sep='\s+', header=0, index_col=0)
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

def addValueToErrorIndicators(s, r, mae, mre, rmse):
    diff = (s - r)
    mae += abs(diff)
    mre += abs(diff / max(1, (s + r)/2))
    rmse += diff**2
    return mae, mre, rmse



def computeErrorRatesFromModelnamesByInterpoliationOptiParalFixedInit(ref, modelname, site_number, startTime, endTime, timestep=1):
    site_name = getSiteNameFromSiteNumber(site_number)
    repoRef = getPathToSimulationDirectoryFromModelname(ref, site_name) 
    print(repoRef)
    repoSimu = getPathToSimulationDirectoryFromModelname(modelname, site_name)


    refHds = fpu.HeadFile(repoRef + '/' + ref + '.hds')
    refTimes = refHds.get_times()
    refKstpkper = refHds.get_kstpkper()
    simuHds = fpu.HeadFile(repoSimu + '/' + modelname + '.hds')
    simuTimes = simuHds.get_times()
    simuKstpkper = simuHds.get_kstpkper()
    

    #Pour chaque jour
    for day in range(startTime, endTime+1):
        print(day)

        # On récupère la matrice de simulation ref
        refHead = refHds.get_data(kstpkper=(0, day))

        nbPeriod = 0
        while (simuTimes[nbPeriod] < day+1) and (nbPeriod < len(simuTimes)):
            nbPeriod+=1
        print("nbPeriod : " + str(nbPeriod))
        #On récupère la matrice de simulation alt supérieure

        print(simuTimes[nbPeriod], day+1)
        if math.isclose(simuTimes[nbPeriod], day+1, rel_tol=1e-3): #simuTimes[nbPeriod] == day+1
            print("condition ==")
            altHeadSup = simuHds.get_data(kstpkper=(timestep-1, nbPeriod))       
            altHeadInf = altHeadSup
            duree = int(simuTimes[nbPeriod])
            pas = day
        else :
            altHeadSup = simuHds.get_data(kstpkper=(timestep-1, nbPeriod))        
            altHeadInf = simuHds.get_data(kstpkper=(timestep-1, nbPeriod-1))
            duree = int(simuTimes[nbPeriod] - simuTimes[nbPeriod-1])
            pas = day - simuTimes[nbPeriod-1]
        

        mae = 0
        mre = 0
        rmse = 0

        for nrow in range(refHead.shape[1]):
            for ncol in range(refHead.shape[2]):
                ss = getNonDryCellHdsValue(altHeadInf, nrow, ncol, refHead.shape[0])
                se = getNonDryCellHdsValue(altHeadSup, nrow, ncol, refHead.shape[0])
                ajoutSimu = (se - ss) / duree
                
                r= getNonDryCellHdsValue(refHead, nrow, ncol,refHead.shape[0])

                s = ss + (ajoutSimu * pas)

                mae, mre, rmse = addValueToErrorIndicators(s, r, mae, mre, rmse)
                
        
        sizeHeads = refHead.shape[1] * refHead.shape[2]
        
        mae = mae / (sizeHeads)
        rmse = math.sqrt(rmse / sizeHeads)
        
        storeErrorValuesIntoCSVFileByInterpolation(ref, modelname, site_name, day, refTimes[day], mae, mre, rmse, startTime, endTime)



def computeFloodsFromModelnamesByInterpoliation(modelname, site_number,timestep=1):
    site_name = getSiteNameFromSiteNumber(site_number)
    repoSimu = getPathToSimulationDirectoryFromModelname(modelname, site_name)
    
    # Get values for topo
    mf = flopy.modflow.Modflow.load(repoSimu + "/" + modelname + '.nam')
    dis = flopy.modflow.ModflowDis.load(repoSimu + "/" + modelname + '.dis', mf)
    topo = dis.top._array
    
    # Get heads values for simulation
    simuHds = fpu.HeadFile(repoSimu + '/' + modelname + '.hds')
    simuTimes = simuHds.get_times()
    #simuKstpkper = simuHds.get_kstpkper()
    
    # Comput parameters
    startTime = 0
    endTime = 15340
    criticalDepth = 0.3

    #endTime+1 - startTime
    floods = {}
    #Pour chaque jour
    for day in range(startTime, endTime+1):
        print(day)

        nbPeriod = 0
        while (simuTimes[nbPeriod] < day+1) and (nbPeriod < len(simuTimes)):
            nbPeriod+=1
        print("nbPeriod : " + str(nbPeriod))
        #On récupère la matrice de simulation alt supérieure

        print(simuTimes[nbPeriod], day+1)
        if math.isclose(simuTimes[nbPeriod], day+1, rel_tol=1e-3): #simuTimes[nbPeriod] == day+1
            print("condition ==")
            altHeadSup = simuHds.get_data(kstpkper=(timestep-1, nbPeriod))       
            altHeadInf = altHeadSup
            duree = int(simuTimes[nbPeriod])
            pas = day
        else :
            altHeadSup = simuHds.get_data(kstpkper=(timestep-1, nbPeriod))        
            altHeadInf = simuHds.get_data(kstpkper=(timestep-1, nbPeriod-1))
            duree = int(simuTimes[nbPeriod] - simuTimes[nbPeriod-1])
            pas = day - simuTimes[nbPeriod-1]

        nbrowtot = altHeadInf.shape[1]        
        nbcoltot = altHeadInf.shape[2]

        onecol = [0]*nbcoltot
        flood = [onecol]
        for i in range(1, nbrowtot):
            flood.append(onecol)

        for nrow in range(nbrowtot):
            for ncol in range(nbcoltot):
                ss = getNonDryCellHdsValue(altHeadInf, nrow, ncol, altHeadInf.shape[0])
                se = getNonDryCellHdsValue(altHeadSup, nrow, ncol, altHeadInf.shape[0])
                ajoutSimu = (se - ss) / duree
                s = ss + (ajoutSimu * pas)

                d = topo[nrow][ncol]- s
                if d > criticalDepth:
                    flood[nrow][ncol]=1
        
        floods[day] = flood
        
    pickle.dump(floods, repoSimu + "/" + modelname + "_floods_pickle")
        
    return floods




def storeErrorValuesIntoCSVFile(ref, modelname, site_name, periodNumber, simulatedDuration, mea, mre, rmse):
    simuRepo = getPathToSimulationDirectoryFromModelname(modelname, site_name)
    if (periodNumber == 0):
        with open(simuRepo + "/" + modelname + '_Ref_' + ref + '_errorsresult.csv', 'w') as f:
            writer = csv.writer(f, delimiter=';')
            writer.writerow(['Period Number', 'Simulated Time', 'MAE', 'MRE', 'RMSE'])
            writer.writerow([periodNumber, simulatedDuration, mea, mre, rmse])
        print("MEA value : ", mea)
        print("MRE value : ", mre)
        print("RMSE value : ", rmse) 


    else:
        with open(simuRepo + "/" + modelname + '_Ref_' + ref + '_errorsresult.csv', 'a') as f:
            writer = csv.writer(f, delimiter=';')
            writer.writerow([periodNumber, simulatedDuration, mea, mre, rmse])
        print("-------------------------")
        print("MEA value : ", mea)
        print("MRE value : ", mre)
        print("RMSE value : ", rmse)

def storeErrorValuesIntoCSVFileByInterpolation(ref, modelname, site_name, periodNumber, simulatedDuration, mea, mre, rmse, startTime, endTime):
    simuRepo = getPathToSimulationDirectoryFromModelname(modelname, site_name)
    nbPart = str(startTime) + "_" + str(endTime)
    if (periodNumber == startTime):
        with open(simuRepo + "/" + modelname + '_Ref_' + ref + '_errorsresult_interpolation_' + str(nbPart) + '.csv', 'w') as f:
            writer = csv.writer(f, delimiter=';')
            writer.writerow(['Period Number', 'Simulated Time', 'MAE', 'MRE', 'RMSE'])
            writer.writerow([periodNumber, simulatedDuration, mea, mre, rmse])
        print("MEA value : ", mea)
        print("MRE value : ", mre)
        print("RMSE value : ", rmse) 


    else:
        with open(simuRepo + "/" + modelname + '_Ref_' + ref + '_errorsresult_interpolation_' + str(nbPart) + '.csv', 'a') as f:
            writer = csv.writer(f, delimiter=';')
            writer.writerow([periodNumber, simulatedDuration, mea, mre, rmse])
        print("-------------------------")
        print("MEA value : ", mea)
        print("MRE value : ", mre)
        print("RMSE value : ", rmse)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-ref", "--ref", type=str, required=False)
    parser.add_argument("-m", "--modelname", type=str, required=True)
    parser.add_argument("-site", "--sitenumber", type=int, required=True)
    parser.add_argument("-ts", "--timestep", type=int, required=False)
    parser.add_argument("-stime", "--starttime", type=int, required=False)
    parser.add_argument("-etime", "--endtime", type=int, required=False)
    parser.add_argument("-flood", "--flood", action='store_true')
    #parser.add_argument("-c", "--corenumber", type=int, required=True)

    args = parser.parse_args()

    ref = args.ref
    modelname = args.modelname
    site_number = args.sitenumber
    timestep = args.timestep
    startTime = args.starttime
    endTime = args.endtime
    flood = args.flood
    #coeur = args.corenumber
    
    if flood:
        computeFloodsFromModelnamesByInterpoliation(modelname, site_number,timestep=1)
    else :
        if (ref is None):
            ref = "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1"

        if timestep is None :
            computeErrorRatesFromModelnamesByInterpoliationOptiParalFixedInit(ref, modelname, site_number, startTime, endTime)
        else :
            computeErrorRatesFromModelnamesByInterpoliationOptiParalFixedInit(ref, modelname, site_number, startTime, endTime, timestep)
    # if timestep is None:
    #     timestep = 1


    # compt=0
    # #coeur=6
    # sims = list(range(1, 15340, 100))
    # sims.append(15340)

    # for ind in range(1, len(sims)):
    #     compt += 1
    #     t = threading.Thread(target=computeErrorRatesFromModelnamesByInterpoliationOptiParal, args=(ref, modelname, site_number, sims[ind-1], sims[ind], timestep))
    #     t.start()
    #     if int(compt / coeur) == compt / coeur:  # Si compt est multiple de 3
    #         t.join()  # alors on attend que les modèles soient terminées pour recommencer
    #         print(compt)
    # t.join() # On attend que les modèles soient finis pour terminer le calcul