
import os
import pandas as pd
import numpy as np
import flopy.modflow as fpm
import flopy.utils as fpu
import numpy.testing as npt
import vtki
import glob
import re
import argparse
import flopy.utils.binaryfile as bf
import csv
import vtk
import sys


def getMeanDdnValuesForWholeSimuTimeForRef(path, nbTimestep, sp, averageSeaLevel, upperLimitForRiskZone):
    meansSum = 0
    dd = getValuesFromHeadFileForATimeStep(
                    path + ".ddn", 9999, 0, 'drawdown')
    dh = getValuesFromHeadFileForATimeStep(path + ".hds", 9999, 0, 'head')
    meansSum += getMeanDdnValuesForATimestepRiskZone(
                    dd, dh, averageSeaLevel, upperLimitForRiskZone)
    for p in range(1, sp) :
        dd = getValuesFromHeadFileForATimeStep(
                    path + ".ddn", nbTimestep-1, p, 'drawdown')
        dh = getValuesFromHeadFileForATimeStep(path + ".hds", nbTimestep-1, p, 'head')
        meansSum += getMeanDdnValuesForATimestepRiskZone(
                    dd, dh, averageSeaLevel, upperLimitForRiskZone)
    return meansSum


def getMeanDdnValuesForWholeSimuTimeRiskZone(path, nbTimestep, sp, averageSeaLevel, upperLimitForRiskZone):
    meansSum = []
    for p in range(0, sp) : # variant : get all kstpkper
        dd = getValuesFromHeadFileForATimeStep(
                    path + ".ddn", nbTimestep-1, p, 'drawdown')
        dh = getValuesFromHeadFileForATimeStep(path + ".hds", nbTimestep-1, p, 'head')
        meansSum.append(getMeanDdnValuesForATimestepRiskZone(
                    dd, dh, averageSeaLevel, upperLimitForRiskZone))
    return np.mean(meansSum)

def getMeanDdnValuesForWholeSimuTime(path, nbTimestep, sp):
    meansSum = []
    for p in range(0, sp) :
        dd = getValuesFromHeadFileForATimeStep(
                    path + ".ddn", nbTimestep-1, p, 'drawdown')
        dh = getValuesFromHeadFileForATimeStep(path + ".hds", nbTimestep-1, p, 'head')
        meansSum.append(getMeanDdnValuesForATimestep(
                    dd))
    return np.mean(meansSum)

## check if iterable only on sp (not timestep and sp)
#Pb for init + simulation, not same sp and timestep

def getValuesFromHeadFileForATimeStep(filePath, timestep, sp, text):
    hdobj = bf.HeadFile(filePath, text=text)
    #print(hdobj.get_kstpkper())
    hds_data = hdobj.get_data(kstpkper=(timestep, sp))
    return hds_data


def getMeanDdnValuesForATimestepRiskZone(dd, dh, averageHdsSealevel, upperLimitForRiskZone):
    nrow = len(dd[0])
    ncol = len(dd[0][0])
    valuesForRiskZone = []
    for j in range(0, nrow):
        for k in range(0, ncol):
            if (dh[0, j, k] != averageHdsSealevel) and (dd[0, j, k] < upperLimitForRiskZone):
                valuesForRiskZone.append(dd[0, j, k])
    return np.mean(valuesForRiskZone)

def getMeanDdnValuesForATimestep(dd):
    nrow = len(dd[0])
    ncol = len(dd[0][0])
    values = []
    for j in range(0, nrow):
        for k in range(0, ncol):
            values.append(dd[0, j, k])
    return np.mean(values)


def getDdnMeansOfAllTimeResultsForASimuFromWatertables(repo):
    os.chdir(repo)
    taille = len(glob.glob("VTU_WaterTable_*.vtu"))

    # Numpy Array to store the mean of ddn data for every stress period
    means_alongTheTime = np.ones(taille)*np.nan

    # One file stores a watertable (ddn data) of the geographical zone (matrix) for a stress period
    for file in glob.glob("VTU_WaterTable_*.vtu"):

        # Retrieving the number of the stress period from the filename
        m = re.search(r'VTU_WaterTable_\w*_(?P<num>\d+).vtu', file)
        if m is not None:
            ind = int(m.group('num'))

        # Retrieving the data from the file
        # Thanks to the vtki library
        gridcolor = vtki.UnstructuredGrid(repo + file)
        nbPointsColor = gridcolor.GetNumberOfCells()
        pointsColor = gridcolor.GetPoints()

        # Retrieving the ddn data for every cell of the geographical zone
        # With the vtki library, it means retreiving the data from the points
        ddn = []
        for i in range(0, nbPointsColor):
            ptcol = pointsColor.GetPoint(i)
            ddn.append(ptcol[2])
        # Transforming the list into a numpy array
        # It will be easier to calculate the means with this type of data structure
        ddn_np = np.array(ddn)
        # Storing the mean for the stress period number 'ind'
        means_alongTheTime[ind] = ddn_np.mean()
    return means_alongTheTime


def getDdnMeansOfAllTimeResultsForASimuFromWatertablesForRiskZonesWithoutSea(repo, averageSeaLevel, upperLimitForRiskZone):
    os.chdir(repo)
    taille = len(glob.glob("VTU_WaterTable_*.vtu"))

    # Numpy Array to store the mean of ddn data for every stress period
    means_alongTheTime = np.ones(taille)*np.nan

    # One file stores a watertable (ddn data) of the geographical zone (matrix) for a stress period
    for file in glob.glob("VTU_WaterTable_*.vtu"):

        # Retrieving the number of the stress period from the filename
        m = re.search(r'VTU_WaterTable_(?P<num>\d+).vtu', file)
        if m is not None:
            ind = int(m.group('num'))

        # Retrieving the ddn data from the file
        # Thanks to the vtki library
        gridcolor = vtki.UnstructuredGrid(repo + file)
        nbPointsColor = gridcolor.GetNumberOfCells()
        pointsColor = gridcolor.GetPoints()

        # Retrieving the ddn data for every cell of the geographical zone
        # With the vtki library, it means retreiving the data from the points
        ddn = []
        for i in range(0, nbPointsColor):
            ptcol = pointsColor.GetPoint(i)
            # Storing the ddn values only if they are below 1m and not equal to the default average sea level
            if (ptcol[2] < upperLimitForRiskZone) and (ptcol[2] != averageSeaLevel):
                ddn.append(ptcol[2])

        # Transforming the list into a numpy array
        # It will be easier to calculate the means with this type of data structure
        ddn_np = np.array(ddn)

        if (ddn_np.size == 0):
            means_alongTheTime[ind] = np.NaN
        else:
            means_alongTheTime[ind] = ddn_np.mean()

    return means_alongTheTime


def getDdnMeansOfAllTimeResultsForASimuFromWatertablesForRiskZones(repo, averageSeaLevel):
    os.chdir(repo)
    taille = len(glob.glob("VTU_WaterTable_*.vtu"))

    # Numpy Array to store the mean of ddn data for every stress period
    means_alongTheTime = np.ones(taille)*np.nan

    # One file stores a watertable (ddn data) of the geographical zone (matrix) for a stress period
    for file in glob.glob("VTU_WaterTable_*.vtu"):

        # Retrieving the number of the stress period from the filename
        m = re.search(r'VTU_WaterTable_(?P<num>\d+).vtu', file)
        if m is not None:
            ind = int(m.group('num'))

        # Retrieving the ddn data from the file
        # Thanks to the vtki library
        gridcolor = vtki.UnstructuredGrid(repo + file)
        nbPointsColor = gridcolor.GetNumberOfCells()
        pointsColor = gridcolor.GetPoints()

        # Retrieving the ddn data for every cell of the geographical zone
        # With the vtki library, it means retreiving the data from the points
        ddn = []
        for i in range(0, nbPointsColor):
            ptcol = pointsColor.GetPoint(i)
            # Storing the ddn values only if they are below 1m and not equal to the default average sea level
            if (ptcol[2] < 1):
                ddn.append(ptcol[2])

        # Transforming the list into a numpy array
        # It will be easier to calculate the means with this type of data structure
        ddn_np = np.array(ddn)

        if (ddn_np.size == 0):
            means_alongTheTime[ind] = np.NaN
        else:
            means_alongTheTime[ind] = ddn_np.mean()

    return means_alongTheTime


def getDdnMeanForWholeSimulationForRiskZonesFromVtuFiles(repo, averageSeaLevel, upperLimitForRiskZone):
    os.chdir(repo)
    taille = len(glob.glob("VTU_WaterTable_*.vtu"))
    
    if (taille == 0):
        print("There is no VTU file in the repository of the given modelname. Check if the files have indeed been generated in the previous phase of the pipeline.")
        sys.exit(0)
    
    # Numpy Array to store the mean of ddn data for every stress period
    means_alongTheTime = np.ones(taille)*np.nan
    # One file stores a watertable (ddn data) of the geographical zone (matrix) for a stress period
    for file in glob.glob("VTU_WaterTable_*.vtu"):
        # Retrieving the number of the stress period from the filename
        m = re.search(r'VTU_WaterTable_*\S*_(?P<num>\d+).vtu', file)
        if m is not None:
            ind = int(m.group('num'))
        else:
            print("There is no VTU file corresponding to the template name in the repository of the given modelname. Check if the files have indeed been generated in the previous phase of the pipeline.")

        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(repo + "/" + file)
        reader.Update() 
        #output = reader.GetOutput()
        heads = reader.GetOutput().GetCellData().GetArray("Heads")
        drawdown = reader.GetOutput().GetCellData().GetArray("Drawdown")

        # Retrieving the ddn data for every cell of the geographical zone
        # With the vtki library, it means retreiving the data from the points
        ddns = []
        for i in range(0, drawdown.GetNumberOfTuples()):
            ddn_val=drawdown.GetTuple(i)[0]
            try:
                npt.assert_approx_equal(heads.GetTuple(i)[0],averageSeaLevel)
                #print("sea zone")
            except:
                if (ddn_val < upperLimitForRiskZone): #(heads.GetTuple(i)[0] !=averageSeaLevel) and 
                    ddns.append(ddn_val)
                elif (heads.GetTuple(i)[0] == averageSeaLevel):
                    print("sea zone 2")
                
        # Transforming the list into a numpy array
        # It will be easier to calculate the means with this type of data structure
        ddn_np = np.array(ddns)
        # Storing the mean for the stress period number 'ind'
        if (ddn_np.size == 0):
            means_alongTheTime[ind] = np.NaN
        else:
            means_alongTheTime[ind] = ddn_np.mean()
    
    return means_alongTheTime

def getDdnMeanForWholeSimulationFromVtuFiles(repo):
    os.chdir(repo)
    taille = len(glob.glob("VTU_WaterTable_*.vtu"))
    if (taille == 0):
        print("There is no VTU file in the repository of the given modelname. Check if the files have indeed been generated in the previous phase of the pipeline.")
        sys.exit(0)
    # Numpy Array to store the mean of ddn data for every stress period
    means_alongTheTime = np.ones(taille)*np.nan

    # One file stores a watertable (ddn data) of the geographical zone (matrix) for a stress period
    for file in glob.glob("VTU_WaterTable_*.vtu"):

        # Retrieving the number of the stress period from the filename
        m = re.search(r'VTU_WaterTable_*\S*_(?P<num>\d+).vtu', file)
        if m is not None:
            ind = int(m.group('num'))

        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(repo + "/" + file)
        reader.Update() 
        drawdown = reader.GetOutput().GetCellData().GetArray("Drawdown")

        # Retrieving the ddn data for every cell of the geographical zone
        # With the vtki library, it means retreiving the data from the points
        ddns = []
        for i in range(0, drawdown.GetNumberOfTuples()):
            ddns.append(drawdown.GetTuple(i)[0])
        # Transforming the list into a numpy array
        # It will be easier to calculate the means with this type of data structure
        ddn_np = np.array(ddns)
        # Storing the mean for the stress period number 'ind'
        if (ddn_np.size == 0):
            means_alongTheTime[ind] = np.NaN
        else:
            means_alongTheTime[ind] = ddn_np.mean()
    
    return means_alongTheTime


def storeMeansValuesIntoCSVFile(dir, modelname, averageSeaLevel, upperLimitForRiskZone, mean, meanall):
    with open(dir + "/" + modelname + '_SeaLvl'+ str(averageSeaLevel) + '_RiskZone' + str(upperLimitForRiskZone) + '_meanresult.csv', 'w') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(['Upper Limit of Ddn values for Risk Zone','level of sea zone removed','Mean', 'Mean all'])
        writer.writerow([upperLimitForRiskZone, averageSeaLevel, mean, meanall])
    print("Mean of ddn values for risk zone ( level <" + str(upperLimitForRiskZone) + ") and without sea region (level = " + str(averageSeaLevel) + "m) : ", mean)
    print("Mean of all ddn values  : ", meanall)

def extractMeansValuesAndStoreThemInCSVFile(modelname, averageSeaLevel, upperLimitForRiskZone):
    #Setting the path to the reprository of the simulation of interest
    repo = "output/simu"
    dir = os.path.join('/'.join(os.getcwd().split('/')[:-1]), repo, modelname)
    # Calculation of drawdown means
    mean = getDdnMeanForWholeSimulationForRiskZonesFromVtuFiles(dir, averageSeaLevel, upperLimitForRiskZone).mean()
    meanall = getDdnMeanForWholeSimulationFromVtuFiles(dir).mean()
    
    # Storing the values into a file
    storeMeansValuesIntoCSVFile(dir, modelname, upperLimitForRiskZone, averageSeaLevel, mean, meanall)


if __name__ == '__main__':
    # Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-sea", "--averagesealevel", type=float, required=False,
                        help="default 0.63m")
    parser.add_argument("-risk", "--riskzonelvl", type=float, required=False,
                        help="default 1m")
    parser.add_argument("-m", "--modelname", type=str, required=True)
    parser.add_argument("-nb", "--nbTimestep", help="value of timestep", type=int, required=False)   
    parser.add_argument("-sp", "--nbStressPeriod", help="number of stress period", type=int, required=False)      
    args = parser.parse_args()

    averageSeaLevel = args.averagesealevel
    upperLimitForRiskZone = args.riskzonelvl
    modelname = args.modelname
    nbTimestep = args.nbTimestep
    sp = args.nbStressPeriod

    extractMeansValuesAndStoreThemInCSVFile(modelname, averageSeaLevel, upperLimitForRiskZone)