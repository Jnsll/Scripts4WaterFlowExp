
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

# Custom imports
import utils



def computeErrorFromGlobalFloodDurationVTUFiles(ref, modelname, upperLimitForFloodZone):
    
    # Get Flood Duration data Matrix from VTU File for Ref
    repoRef = utils.getPathToSimulationDirectoryFromModelname(ref)
    #os.chdir(repoRef)

    refVTU = utils.getFloodDurationVTUFileNameFromModelnameAndLimitValueForFloodZone(ref, upperLimitForFloodZone)
    ## Need to check if file exists ! Else : sys.exit(0)
    print(refVTU)
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(repoRef + "/" + refVTU)
    reader.Update() 
    floodDurationRefData = reader.GetOutput().GetCellData().GetArray("FloodDuration")
    tailleRef = floodDurationRefData.GetNumberOfTuples()

    # for i in range(0, floodDurationRefData.GetNumberOfTuples()):
    #     flood_valRef = floodDurationRefData.GetTuple(i)[0]
    repoExpe = utils.getPathToSimulationDirectoryFromModelname(modelname)
    expeVTU = utils.getFloodDurationVTUFileNameFromModelnameAndLimitValueForFloodZone(modelname, upperLimitForFloodZone)
    print(expeVTU)
    reader.SetFileName(repoExpe + "/" + expeVTU)
    reader.Update()
    floodDurationExpeData = reader.GetOutput().GetCellData().GetArray("FloodDuration")
    tailleExpe = floodDurationExpeData.GetNumberOfTuples()

    print("ref: " + str(tailleRef))
    print("expe: " + str(tailleExpe))

    mae = 0
    mre = 0
    rmse = 0

    for i in range(tailleRef):
        diff = floodDurationRefData.GetTuple(i)[0] - floodDurationExpeData.GetTuple(i)[0]
        mae += abs(diff)
        mre += abs(diff / max(1, (floodDurationRefData.GetTuple(i)[0] + floodDurationExpeData.GetTuple(i)[0])/2))
        rmse += diff**2

    mae = mae / tailleRef
    rmse = math.sqrt(rmse / tailleRef)

    return mae, mre, rmse
    # print("mea : " + str(mae))
    # print("mre : " + str(mre))   
    # print("rmse : " + str(rmse))


def storeErrorValuesIntoCSVFile(ref, modelname, averageSeaLevel, upperLimitForFloodZone, mea, mre, rmse):
    repoExpe = utils.getPathToSimulationDirectoryFromModelname(modelname)
    with open(repoExpe + "/" + modelname + '_Ref_' + ref + '_SeaLvl'+ str(averageSeaLevel) + '_FloodZone' + str(upperLimitForFloodZone) + '_errorsresult.csv', 'w') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(['Upper Limit of thickness values for Flood Zone','Level of removed sea zone','MAE', 'MRE', 'RMSE'])
        writer.writerow([upperLimitForFloodZone, averageSeaLevel, mea, mre, rmse])
    print("MEA value for flood zone ( level <" + str(upperLimitForFloodZone) + "m) and without sea region (level = " + str(averageSeaLevel) + "m) : ", mea)
    print("MRE value for flood zone ( level <" + str(upperLimitForFloodZone) + "m) and without sea region (level = " + str(averageSeaLevel) + "m) : ", mre)
    print("RMSE value for flood zone ( level <" + str(upperLimitForFloodZone) + "m) and without sea region (level = " + str(averageSeaLevel) + "m) : ", rmse)


def computeErrorsAndStoreThemIntoCSVFile(ref, modelname, averageSeaLevel, upperLimitForFloodZone):
    mae, mre, rmse = computeErrorFromGlobalFloodDurationVTUFiles(ref, modelname, upperLimitForFloodZone)
    storeErrorValuesIntoCSVFile(ref, modelname, averageSeaLevel, upperLimitForFloodZone, mae, mre, rmse)


if __name__ == '__main__':
    # Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-sea", "--averagesealevel", type=float, required=False,
                        help="default 0.63m")
    parser.add_argument("-flood", "--floodzonelvl", type=float, required=False,
                        help="default 0.1m")
    parser.add_argument("-m", "--modelname", type=str, required=True)
    parser.add_argument("-ref", "--reference", type=str, required=True)
    
   
    args = parser.parse_args()

    averageSeaLevel = args.averagesealevel
    upperLimitForFloodZone = args.floodzonelvl
    modelname = args.modelname
    ref = args.reference
    
    computeErrorsAndStoreThemIntoCSVFile(ref, modelname, averageSeaLevel, upperLimitForFloodZone)

    


