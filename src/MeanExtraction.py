
import os
import pandas as pd
import numpy as np
import flopy.modflow as fpm
import flopy.utils as fpu


def getDataFromSimuResultFile(direct, modelname, extension):
    hfile = fpu.HeadFile(os.path.join(direct, modelname + extension), text='drawdown')
    h = hfile.get_data()
    return h


# directory = "/DATA/These/Projects/WaterFlowSimulationModel/runWithInitStressPeriod/"
# modelname= "model_base_multi_stressPeriod"
# extension = ".ddn"
# data_ddn = getDataFromSimuResultFile(directory, modelname, extension)
# df_ddn = pd.DataFrame(data=np.array(data_ddn[0]))
# print(df_ddn.shape)
# print(df_ddn.describe())
#df_mean_ddn = pd.DataFrame(df_ddn.mean())





import vtki
import glob, os
import re
import argparse

### Parser ###
parser = argparse.ArgumentParser()
parser.add_argument("dir", help="Path of the directory containing all the watertable files", type=str)
args = parser.parse_args()
repo = args.dir #"/DATA/These/Projects/WaterFlowSimulationModel/runWithInitStressPeriod/"

def getDdnMeansOfAllTimeResultsForASimuFromWatertables(repo):
    os.chdir(repo)
    taille=len(glob.glob("VTU_WaterTable_*.vtu"))
    means_alongTheTime = np.empty(taille)

    for file in glob.glob("VTU_WaterTable_*.vtu"):
        #print(file)
        m = re.search(r'VTU_WaterTable_(?P<num>\d+).vtu', file)
        if m is not None:
            ind= int(m.group('num'))
        filenameColor = file
        gridcolor = vtki.UnstructuredGrid(repo + filenameColor)
        nbPointsColor = gridcolor.GetNumberOfCells()
        pointsColor = gridcolor.GetPoints()
    #print(pointsColor)
        ddn = []
        for i in range(0,nbPointsColor) :
            ptcol = pointsColor.GetPoint(i)
            ddn.append(ptcol[2])
    #print(ptcol[2])

        ddn_np = np.array(ddn) 
        means_alongTheTime[ind] = ddn_np.mean()
    return means_alongTheTime


if __name__ == '__main__':
    means = getDdnMeansOfAllTimeResultsForASimuFromWatertables(repo)
    print(means.mean())