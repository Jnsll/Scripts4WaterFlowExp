
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




def insertFloodZoneValuesForEachPeriodIntoVtuFiles(repo, modelname, averageSeaLevel, upperLimitForFloodZone, refSimu):
    os.chdir(repo)
    vtuFiles = glob.glob("VTU_WaterTable_*[!.][0-9].vtu")

    if (len(vtuFiles) == 0):
        print("There is no VTU file in the repository of the given modelname. Check if the files have indeed been generated in the previous phase of the pipeline.")
        sys.exit(0)
    
    # One file stores a watertable of the geographical zone (matrix) for a stress period
    
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(repo + "/" + vtuFiles[0])
    reader.Update() 
    heads = reader.GetOutput().GetCellData().GetArray("Heads")
    taille = heads.GetNumberOfTuples()
    floodDuration = np.zeros(taille)

    #shutil.copyfile(vtuFiles[0], "VTU_WaterTable_FloodDuration.vtu")
    fin = open(vtuFiles[0], "r" )
    data_list = fin.readlines()
    fin.close()

    del data_list[4:3522]

    fout = open("VTU_WaterTable_" + modelname + "_FloodDuration_" + str(upperLimitForFloodZone) + ".vtu", "w")
    fout.writelines(data_list)
    fout.close()

    hds = fpu.HeadFile(modelname + '.hds')
    times = hds.get_times()

    for file in vtuFiles:
        m = re.search(r'VTU_WaterTable_*\S*_(?P<num>\d+).vtu', file)
        if m is not None:
            ind = int(m.group('num'))
        else:
            print("There is no VTU file corresponding to the template name in the repository of the given modelname. Check if the files have indeed been generated in the previous phase of the pipeline.")


        print(file)
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(repo + "/" + file)
        reader.Update() 
        #output = reader.GetOutput()
        flood = reader.GetOutput().GetCellData().GetArray("FloodZone")
        # if flood is not None:
        #     print("Flood zones have been already inserted into the vtu file named : " + file)

        heads = reader.GetOutput().GetCellData().GetArray("Heads")
        drawdown = reader.GetOutput().GetCellData().GetArray("Drawdown")
        
        tailleb = drawdown.GetNumberOfTuples()
        # Retrieving the ddn data for every cell of the geographical zone
        # With the vtki library, it means retreiving the data from the points
        ddns = np.zeros(tailleb)
        for i in range(0, drawdown.GetNumberOfTuples()):
            ddn_val=drawdown.GetTuple(i)[0]
            try:
                npt.assert_approx_equal(heads.GetTuple(i)[0],averageSeaLevel)
                #print("sea zone")
            except:
                if (ddn_val < upperLimitForFloodZone): #(heads.GetTuple(i)[0] !=averageSeaLevel) and 
                        if (ind == 0):
                            if refSimu is None:
                                ddns[i] = round(times[ind])
                        else :
                            ddns[i] = round(times[ind] - times[ind-1])
                elif (heads.GetTuple(i)[0] == averageSeaLevel):
                    print("sea zone 2")
                
        # Transforming the list into a numpy array
        # It will be easier to calculate the means with this type of data structure
        ddn_np = np.array(ddns)
        # Storing the mean for the stress period number 'ind'
        #f = File(repo + model + "_" + ind + ".csv")
        if flood is None :

            f = open(file, "r")
            contents = f.readlines()
            f.close()

            value= '        <DataArray type="Float64" Name="FloodZone" format="ascii">\n'
            for item in range(len(ddn_np)):
                textvalue = str(ddn_np[item])
                if item == 0:
                    value += '          ' + textvalue + ' '
                elif item % 20 == 0:
                    value+= textvalue + '\n          '
                else:
                    value+= textvalue + ' '
            value+= '\n'
            value+= '        </DataArray>\n'
            contents.insert(4, value)

            f = open(file, "w")
            contents = "".join(contents)
            f.write(contents)
            f.close()
        else:
            print("Flood zones have been already inserted into the vtu file named : " + file)
            
            l = [floodDuration, ddn_np]
            floodDuration = np.sum(l, axis=0)

    return floodDuration

def extractFloodZoneDurationFromVTUFiles(modelname, averageSeaLevel, upperLimitForFloodZone, refSimu):

    repo = "output/simu"
    dir = os.path.join('/'.join(os.getcwd().split('/')[:-1]), repo, modelname)

    floodDuration = insertFloodZoneValuesForEachPeriodIntoVtuFiles(dir, modelname, averageSeaLevel, upperLimitForFloodZone, refSimu)
    
    f = open("VTU_WaterTable_" + modelname + "_FloodDuration_" + str(upperLimitForFloodZone) + ".vtu", "r")
    contents = f.readlines()
    f.close()

    value= '        <DataArray type="Float64" Name="FloodDuration" format="ascii">\n'
    for item in range(len(floodDuration)):
        textvalue = str(floodDuration[item])
        if item == 0:
            value += '          ' + textvalue + ' '
        elif item % 20 == 0:
            value+= textvalue + '\n          '
        else:
            value+= textvalue + ' '
    value+= '\n'
    value+= '        </DataArray>\n'
    contents.insert(4, value)

    f = open("VTU_WaterTable_" + modelname + "_FloodDuration_" + str(upperLimitForFloodZone) + ".vtu", "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()



if __name__ == '__main__':
    # Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-sea", "--averagesealevel", type=float, required=False,
                        help="default 0.63m")
    parser.add_argument("-flood", "--floodzonelvl", type=float, required=False,
                        help="default 0.1m")
    parser.add_argument("-m", "--modelname", type=str, required=True)
    parser.add_argument("-r", "--refsimu", action='store_true')
   
    args = parser.parse_args()

    averageSeaLevel = args.averagesealevel
    upperLimitForFloodZone = args.floodzonelvl
    modelname = args.modelname
    refSimu = args.refsimu
 
    extractFloodZoneDurationFromVTUFiles(modelname, averageSeaLevel, upperLimitForFloodZone, refSimu)

    


