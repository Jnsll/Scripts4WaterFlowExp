# -*- coding: utf-8 -*-
"""
    Provides generic functions used in several scripts.
"""

import os
import flopy
import pandas as pd

__author__ = "June Sallou"
__maintainer__ = "June Sallou"
__credits__ = ["June Sallou"]
__license__ = "MIT"
__version__ = "0.0.3"
__date__ = "04/30/2019"
__email__ = "june.benvegnu-sallou@univ-rennes1.fr"



mainAppRepo = os.path.dirname(os.path.abspath(__file__)) + '/' # A verifier , il faudrait surement remonter d'un rang superieur

def get_model_name(site_number, chronicle, approx, rate, ref):
    model_name = "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_site" + str(site_number) + "_Chronicle" + str(chronicle)
    if not ref:
        model_name += "_Approx" + str(approx)
        if approx == 0:
            model_name += "_Period" + str(rate)
        elif approx==1:
            model_name += "_RechThreshold" + str(rate)
    return model_name

def get_site_name_from_site_number(site_number):
    sites = pd.read_csv(mainAppRepo + 'data/study_sites.txt',
                        sep=',', header=0, index_col=0) #\\s+
    site_name = sites.index._data[site_number] + '/'
    return site_name

def get_path_to_simulation_directory(site_number, chronicle, approx, rate, ref):
    model_name = get_model_name(site_number, chronicle, approx, rate, ref)
    site_name = get_site_name_from_site_number(site_number)
    return os.path.join(mainAppRepo, "outputs/", site_name, model_name)


def importDataFromModel(modelname, dataToLoad):
    """ 
        param dataToLoad : a list with the type of data to load from the model
        example : ['upw', 'dis']
    """
    mf1 = flopy.modflow.Modflow.load(modelname + '.nam', verbose=False,check=False, load_only=dataToLoad)
    return mf1


def writeExecutionTimeInLogfile(path, modelname, duration):
    with open(path + '/' + modelname + '_log.txt', 'w') as f: 
        f.write('Execution time (s) of ' + modelname +'\n')
        f.write(str(duration))

def getPathToSimulationDirectoryFromModelname(modelname):
    repo = "data"
    return os.path.join('/'.join(os.path.realpath(__file__).split('/')[:-2]), repo, modelname) 

def getFloodDurationVTUFileNameFromModelnameAndLimitValueForFloodZone(modelname, upperLimitForFloodZone):
    return "VTU_WaterTable_" + modelname + "_FloodDuration_" + str(upperLimitForFloodZone) + ".vtu"


def generate_model_name(chronicle, approx, rate, ref, steady, site=None, time_param=0, geology_param=0, thickness_param=1, permeability_param=86.4, theta_param=0.1, step=None):
        model_name = r'model_' + 'time_' + str(time_param) + '_geo_' + str(geology_param) + '_thick_' + str(thickness_param)
        if geology_param == 0:
            model_name = model_name + '_K_' + str(permeability_param) + '_Sy_' + str(theta_param)
        if step is None:
            model_name += "_Step" + str(1)
        else:
            model_name += "_Step" + str(step)
        if site is not None:
            model_name += "_site" + str(site)
        if chronicle is not None:
            model_name += "_Chronicle" + str(chronicle)
        else:
            model_name += "_Chronicle" + str(chronicle)
        if (not ref) and approx is not None:
            model_name += "_Approx" + str(approx)
            if approx==0:
                model_name += "_Period" + str(rate)
            elif approx==1:
                model_name += "_RechThreshold" + str(rate)
        if steady:
            model_name += "_SteadyState"
        return model_name