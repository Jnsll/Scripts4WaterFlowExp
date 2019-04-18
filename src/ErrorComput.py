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


mainAppRepo = os.getcwd() + "/"


def getPathToSimulationDirectoryFromModelname(modelname, site_name):
    return os.path.join(mainAppRepo, site_name, modelname)

def getSiteNameFromSiteNumber(site_number):
    sites = pd.read_csv(mainAppRepo + "data/study_sites.txt", sep='\s+', header=0, index_col=0)
    site_number = 2
    site_name = sites.index._data[site_number] + '/'
    return site_name


def computeErrorRatesFromModelnames(ref, modelname, site_number, timestep=1):
    site_name = getSiteNameFromSiteNumber(site_number)

    #repoRef = "/DATA/These/Publi/eScience19/Expes/model_ref/Agon-Coutainville/model_time_0_geo_0_thick_1_K_86.4_Sy_0.1/"
    repoRef = getPathToSimulationDirectoryFromModelname(ref, site_name) 
    #repoSimu = "/DATA/These/Projects/Model/docker/app/Agon-Coutainville/model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodsemester_Step1"
    repoSimu = getPathToSimulationDirectoryFromModelname(modelname, site_name) 

    refHds = fpu.HeadFile(repoRef + '/' + ref + '.hds')
    refTimes = refHds.get_times()
    refKstpkper = refHds.get_kstpkper()
    simuHds = fpu.HeadFile(repoSimu + '/' + modelname + '.hds')
    simuTimes = simuHds.get_times()
    simuKstpkper = simuHds.get_kstpkper()

    for numPrd in range(0, len(simuTimes)): #On prend l'initialisation
        if (numPrd == 0):
            refHead = refHds.get_data(kstpkper=(0, simuTimes[numPrd]-1))
        else:
            refHead = refHds.get_data(kstpkper=(timestep-1, simuTimes[numPrd]-1))
        simuHead = simuHds.get_data(kstpkper=simuKstpkper[numPrd])

        mae = 0
        mre = 0
        rmse = 0

        for nrow in range(simuHead.shape[1]):
            for ncol in range(simuHead.shape[2]):
                layerS = 0
                layerR = 0
                s = simuHead[0][nrow][ncol]
                r = refHead[0][nrow][ncol]
                while (math.isclose(abs(s)/1e+30, 1, rel_tol=1e-3)):
                    s = simuHead[layerS-1][nrow][ncol]
                    layerS-=1
                while (math.isclose(abs(r)/1e+30, 1, rel_tol=1e-3)):
                    r = refHead[layerR-1][nrow][ncol]
                    layerR-=1
                diff = s - r
                mae += abs(diff)
                mre += abs(diff / max(1, (simuHead[0][nrow][ncol] - refHead[0][nrow][ncol])/2))
                rmse += diff**2
        
        sizeHeads = simuHead.shape[1] * simuHead.shape[2]
        
        mae = mae / (sizeHeads)
        rmse = math.sqrt(rmse / sizeHeads)

        storeErrorValuesIntoCSVFile(ref, modelname, site_name, numPrd, simuTimes[numPrd], mae, mre, rmse)
        




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

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-ref", "--ref", type=str, required=False)
    parser.add_argument("-m", "--modelname", type=str, required=True)
    parser.add_argument("-site", "--sitenumber", type=int, required=True)
    parser.add_argument("-ts", "--timestep", type=int, required=False)
    args = parser.parse_args()

    ref = args.ref
    modelname = args.modelname
    site_number = args.sitenumber
    timestep = args.timestep

    if ref is None:
        ref = "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1"

    if timestep is None :
        computeErrorRatesFromModelnames(ref, modelname, site_number)
    else :
        computeErrorRatesFromModelnames(ref, modelname, site_number, timestep)
    


