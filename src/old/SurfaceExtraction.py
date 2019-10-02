
import os
import pandas as pd
import numpy as np
import flopy.modflow as fpm
import flopy.utils as fpu

import vtki
import glob, os
import re
import argparse

### Parser ###
parser = argparse.ArgumentParser()
parser.add_argument("dir", help="Path of the directory containing all the watertable files", type=str)


args = parser.parse_args()
repo = args.dir #"/DATA/These/Projects/WaterFlowSimulationModel/runWithInitStressPeriod/"
typeOfCalculation = args.type

averageSeaLevel = 0.63


def getDdnMeansOfAllTimeResultsForASimuFromWatertables(repo):
    os.chdir(repo)
    taille=len(glob.glob("VTU_WaterTable_*.vtu"))

    #Numpy Array to store the mean of ddn data for every stress period
    means_alongTheTime = np.ones(taille)*np.nan

    # One file stores a watertable (ddn data) of the geographical zone (matrix) for a stress period
    for file in glob.glob("VTU_WaterTable_*.vtu"):

        #Retrieving the number of the stress period from the filename
        m = re.search(r'VTU_WaterTable_(?P<num>\d+).vtu', file)
        if m is not None:
            ind= int(m.group('num'))

        #Retrieving the data from the file 
        #Thanks to the vtki library
        gridcolor = vtki.UnstructuredGrid(repo + file)
        nbPointsColor = gridcolor.GetNumberOfCells()
        pointsColor = gridcolor.GetPoints()

        # Retrieving the ddn data for every cell of the geographical zone 
        # With the vtki library, it means retreiving the data from the points
        ddn = []
        for i in range(0,nbPointsColor) :
            ptcol = pointsColor.GetPoint(i)
            ddn.append(ptcol[2])
        #Transforming the list into a numpy array
        #It will be easier to calculate the means with this type of data structure
        ddn_np = np.array(ddn) 
        #Storing the mean for the stress period number 'ind'
        means_alongTheTime[ind] = ddn_np.mean()
    return means_alongTheTime

def getDdnMeansOfAllTimeResultsForASimuFromWatertablesForRiskZonesWithoutSea(repo, averageSeaLevel):
    os.chdir(repo)
    taille=len(glob.glob("VTU_WaterTable_*.vtu"))

    #Numpy Array to store the mean of ddn data for every stress period
    means_alongTheTime = np.ones(taille)*np.nan

    # One file stores a watertable (ddn data) of the geographical zone (matrix) for a stress period
    for file in glob.glob("VTU_WaterTable_*.vtu"):
        
        #Retrieving the number of the stress period from the filename
        m = re.search(r'VTU_WaterTable_(?P<num>\d+).vtu', file)
        if m is not None:
            ind= int(m.group('num'))
        
        #Retrieving the ddn data from the file 
        #Thanks to the vtki library
        gridcolor = vtki.UnstructuredGrid(repo + file)
        nbPointsColor = gridcolor.GetNumberOfCells()
        pointsColor = gridcolor.GetPoints()

        # Retrieving the ddn data for every cell of the geographical zone 
        # With the vtki library, it means retreiving the data from the points
        ddn = []
        for i in range(0,nbPointsColor) :
            ptcol = pointsColor.GetPoint(i)
            #Storing the ddn values only if they are below 1m and not equal to the default average sea level
            if (ptcol[2]<1) and (ptcol[2]!=averageSeaLevel):
                ddn.append(ptcol[2])
        
        #Transforming the list into a numpy array
        #It will be easier to calculate the means with this type of data structure
        ddn_np = np.array(ddn) 
        # if (ddn_np.size ==0):
        #     means_alongTheTime[ind] = np.NaN
        # else:
        if (ddn_np.size ==0):
            means_alongTheTime[ind] = np.NaN
        else:
            means_alongTheTime[ind] = ddn_np.mean()        

    return means_alongTheTime

def getDdnMeansOfAllTimeResultsForASimuFromWatertablesForRiskZones(repo, averageSeaLevel):
    os.chdir(repo)
    taille=len(glob.glob("VTU_WaterTable_*.vtu"))

    #Numpy Array to store the mean of ddn data for every stress period
    means_alongTheTime = np.ones(taille)*np.nan

    # One file stores a watertable (ddn data) of the geographical zone (matrix) for a stress period
    for file in glob.glob("VTU_WaterTable_*.vtu"):
        
        #Retrieving the number of the stress period from the filename
        m = re.search(r'VTU_WaterTable_(?P<num>\d+).vtu', file)
        if m is not None:
            ind= int(m.group('num'))
        
        #Retrieving the ddn data from the file 
        #Thanks to the vtki library
        gridcolor = vtki.UnstructuredGrid(repo + file)
        nbPointsColor = gridcolor.GetNumberOfCells()
        pointsColor = gridcolor.GetPoints()

        # Retrieving the ddn data for every cell of the geographical zone 
        # With the vtki library, it means retreiving the data from the points
        ddn = []
        for i in range(0,nbPointsColor) :
            ptcol = pointsColor.GetPoint(i)
            #Storing the ddn values only if they are below 1m and not equal to the default average sea level
            if (ptcol[2]<1):
                ddn.append(ptcol[2])
        
        #Transforming the list into a numpy array
        #It will be easier to calculate the means with this type of data structure
        ddn_np = np.array(ddn) 
        # if (ddn_np.size ==0):
        #     means_alongTheTime[ind] = np.NaN
        # else:
        if (ddn_np.size ==0):
            means_alongTheTime[ind] = np.NaN
        else:
            means_alongTheTime[ind] = ddn_np.mean()        

    return means_alongTheTime



if __name__ == '__main__':
    if typeOfCalculation=="all":
        means = getDdnMeansOfAllTimeResultsForASimuFromWatertables(repo)
        print(means.mean())
    elif typeOfCalculation=="riskSea":
        means = getDdnMeansOfAllTimeResultsForASimuFromWatertablesForRiskZonesWithoutSea(repo, averageSeaLevel)
        if (means[~np.isnan(means)].size != 0):
            print(means[~np.isnan(means)].mean())
        else:
            print("No value")
    elif typeOfCalculation=="riskAll":
        means = getDdnMeansOfAllTimeResultsForASimuFromWatertablesForRiskZones(repo, averageSeaLevel)
        print(means[~np.isnan(means)].mean())
    else :
        print("The type parameter can be either 'all' or 'risk'. Please enter one of them only.")