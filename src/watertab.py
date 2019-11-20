import sys, os
import pandas as pd
import argparse
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'vtk_export_watertable'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'custom_utils'))
from custom_utils import helpers as utils
from vtk_export_watertable import vtk_export_watertable as vtk_watertable





def get_model_name(site_number, chronicle, approx, rate, ref):
    model_name = "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_site" + str(site_number) + "_Chronicle" + str(chronicle)
    if not ref:
        model_name += "_Approx" + str(approx)
        if approx == 0:
            model_name += "_Period" + str(rate)
        elif approx==1:
            model_name += "_RechThreshold" + str(rate)
    return model_name


def watertab(site, chronicle, approx, rate, ref, folder):
    folder_path = os.path.dirname(os.path.abspath(__file__)) + '/'
    sites = pd.read_table(folder_path + "data/study_sites.txt", sep=',', header=0, index_col=0)
    model_name = get_model_name(site, chronicle, approx, rate, ref)
    site_name = sites.index._data[site]
    mainRepo = folder + '/' + site_name + '/' + model_name + '/'
    coordinates = sites._get_values[site,1:5]

    vtk_watertable(modelname=model_name, modelfolder=mainRepo, coord=coordinates)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-site", "--site", type=int, help= "2: Agon-Coutainville or 3:Saint-Germain-Sur-Ay", required=True)
    parser.add_argument("-approx", "--approximation", type=int, required=True)
    parser.add_argument("-chr", "--chronicle", type=int)
    parser.add_argument("-rate", "--rate", type=float, required=False)
    parser.add_argument("-f", "--folder", type=str, required=True)
    parser.add_argument("-ref", "--ref", action='store_true')
    args = parser.parse_args()


    site = args.site
    chronicle = args.chronicle
    folder= args.folder
    approx = args.approximation
    chronicle = args.chronicle
    rate = args.rate
    ref = args.ref

    watertab(site, chronicle, approx, rate, ref, folder)