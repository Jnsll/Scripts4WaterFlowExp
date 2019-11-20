# -*- coding: utf-8 -*-
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


mainAppRepo = os.path.dirname(os.path.abspath(__file__)) + '/'

# Simulation features
nb_years = 42 # Duration
startTime = 0
endTime = 15340 # end of the 42 years of the simulation
# H indicator computation features
tsmax = 121  # 4 months = 121,75 days
max_diff_threshold = 7 # max discrancy of the duration of flood episode
dc = 0.3  # depth of threshold for the undeground vulnerable zone
alpha = float(1/3)


def cut_hds_ref_file_into_one_day_hds_files(site_number, chronicle, approx, rate, folder):
    ref_name = get_model_name(site_number, chronicle, approx, rate, ref=True)
    if folder is None:
        repo_ref = get_path_to_simulation_directory(site_number, chronicle, approx, rate, ref=True)
    else :
        site_name = get_site_name_from_site_number(site_number)
        repo_ref = folder + site_name + ref_name

    ref_hds = fpu.HeadFile(repo_ref + '/' + ref_name + '.hds')  
    os.makedirs(repo_ref + "/HDS", exist_ok=True)

    for day in range(startTime, endTime+1):
        refHead = ref_hds.get_data(kstpkper=(0, day))
        np.save(repo_ref + '/' + 'HDS' + '/' + 'hds_' + str(day) + '.npy', refHead)
        print((day/endTime)*100, "%")

def get_model_name(site_number, chronicle, approx, rate, ref):
    model_name = "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_site" + str(site_number) + "_Chronicle" + str(chronicle)
    if not ref:
        model_name += "_Approx" + str(approx)
        if approx == 0:
            model_name += "_Period" + str(rate)
        elif approx==1:
            model_name += "_RechThreshold" + str(rate)
    return model_name


def get_path_to_simulation_directory(site_number, chronicle, approx, rate, ref):
    model_name = get_model_name(site_number, chronicle, approx, rate, ref)
    site_name = get_site_name_from_site_number(site_number)
    return os.path.join(mainAppRepo, "outputs/", site_name, model_name)


def get_site_name_from_site_number(site_number):
    sites = pd.read_csv(mainAppRepo + 'data/study_sites.txt',
                        sep=',', header=0, index_col=0) #\\s+
    site_name = sites.index._data[site_number] + '/'
    return site_name


def get_non_dry_cell_hds_value(hds, nrow, ncol, nlayer):
    layer = 0
    h = hds[layer][nrow][ncol]
    while (math.isclose(abs(h)/1e+30, 1, rel_tol=1e-3)) and layer < nlayer:
        if layer == nlayer-1:
            print("cell completely dry")
        else:
            h = hds[layer+1][nrow][ncol]
            layer += 1
    return h


def getWeightToSurface(zs, h, dc, alpha):
    """
        zs : value of soil surface at the point x
        h : head value for watertable at the x point
        dc : critical depth
        alpha : ratio of critical depth for calculating the width of transition zone
    """
    ddc = alpha * dc  # Alpha must be not null

    borneInf = zs - (dc + (ddc / 2))
    borneSup = zs - (dc - (ddc / 2))

    if h <= borneInf:
        Ws = 0
    elif h >= borneSup:
        Ws = 1
    else:
        Ws = math.sin((math.pi * (h - borneSup)) / (2*ddc))

    return Ws


def getNonDryCellHdsValue(hds, nrow, ncol, nlayer):
    layer = 0
    h = hds[layer][nrow][ncol]
    while (math.isclose(abs(h)/1e+30, 1, rel_tol=1e-3)) and layer < nlayer:
        if layer == nlayer-1:
            print("cell completely dry")
        else:
            h = hds[layer+1][nrow][ncol]
            layer += 1
    return h

def get_soil_surface_values_for_a_simulation(repo_simu, model_name):
    """
        Retrieve the matrix of the topographical altitude.
    """
    mf = flopy.modflow.Modflow.load(repo_simu + "/" + model_name + '.nam')
    dis = flopy.modflow.ModflowDis.load(
        repo_simu + "/" + model_name + '.dis', mf)
    topo = dis.top._array
    
    return topo

def topo_file(site_number, chronicle, approx, rate, ref, folder):
    model_name = get_model_name(site_number, chronicle, approx, rate, ref)
    site_name = get_site_name_from_site_number(site_number)
    repo = folder + site_name + model_name
    topo = get_soil_surface_values_for_a_simulation(repo, model_name)
    np.save(repo + "/soil_surface_topo_"+ model_name + ".npy", topo)
    print(repo + "/soil_surface_topo_"+ model_name + ".npy")


def compute_h_error_by_interpolation(site_number, chronicle, approx, rate, ref, folder, time_step=1):
    # Get data for reference simulation
    # Path to repository
    ref_name = get_model_name(site_number, chronicle, approx, rate, ref=True)
    if folder is None:
        repo_ref = get_path_to_simulation_directory(site_number, chronicle, approx, rate, ref=True)
    else :
        site_name = get_site_name_from_site_number(site_number)
        repo_ref = folder + site_name + ref_name

    topo_ref = np.load(repo_ref + "/soil_surface_topo_"+ ref_name + ".npy")        
    # Watertable altitude
    ref_hds = fpu.HeadFile(repo_ref + '/' + ref_name + '.hds')    
    
    if ref:
        # Get data for alternate simulation
        # Path to repository
        repo_simu = repo_ref 
        simu_name = ref_name
        # Topographical altitude
        topoSimu = topo_ref
        # Get heads values for simulation
        simu_hds = ref_hds
        
    else:
        # Get data for alternate simulation
        # Path to repository
        simu_name = get_model_name(site_number, chronicle, approx, rate, ref=False)
        if folder is None:
            repo_simu = get_path_to_simulation_directory(site_number, chronicle, approx, rate, ref=False)
            
        else:
            repo_simu = folder + site_name + simu_name
        # Topographical altitude
        topoSimu = get_soil_surface_values_for_a_simulation(repo_simu, simu_name)
        # Get heads values for simulation
        simu_hds = fpu.HeadFile(repo_simu + '/' + simu_name + '.hds')

    # Duration of simulated periods
    simu_times = simu_hds.get_times()

    # Initialisation of error computing variables
    smWs = 0
    sherrorsup = 0

    # For every day, we have to compare the reference and alternate simulations
    for day in range(startTime, endTime+1):
        print((day/endTime)*100, "%")

        # Retrieve the watertable altitude values for the reference simulation for the considered day
        refHead = ref_hds.get_data(kstpkper=(0, day))

        # Compute the number of the simulated period of the alternate simulation for the considered day
        nbPeriod = 0
        while (simu_times[nbPeriod] < day+1) and (nbPeriod < len(simu_times)):
            nbPeriod += 1
        # print("nbPeriod : " + str(nbPeriod))

        # When the considered day match with the duration of the corresponding simulated period
        # We retrieve the value of the watertable altitude from the corresponding matrix
        if math.isclose(simu_times[nbPeriod], day+1, rel_tol=1e-3):
            altHeadSup = simu_hds.get_data(kstpkper=(time_step-1, nbPeriod))
            altHeadInf = altHeadSup
            duree = int(simu_times[nbPeriod])
            pas = 0
        # Otherwise, we have to interpolate the watertable altitude value for the considered day
        else:
            # The considered day is situated between the simulated period number 'nbPeriod-1' and number 'nbPeriod'
            altHeadSup = simu_hds.get_data(kstpkper=(time_step-1, nbPeriod))
            altHeadInf = simu_hds.get_data(kstpkper=(time_step-1, nbPeriod-1))
            duree = int(simu_times[nbPeriod] - simu_times[nbPeriod-1])
            pas = day - simu_times[nbPeriod-1]

        nbrowtot = altHeadInf.shape[1]
        nbcoltot = altHeadInf.shape[2]

        # We want to store the presence of cells part of a flood episode
        # We go through all the cells of the matrix representing the study site
        for nrow in range(nbrowtot):
            for ncol in range(nbcoltot):
                # Watertable altitude value for simulated period with duration lower than considered day
                ss = get_non_dry_cell_hds_value(
                    altHeadInf, nrow, ncol, altHeadInf.shape[0])
                # Watertable altitude value for simulated period with duration higher than considered day
                se = get_non_dry_cell_hds_value(
                    altHeadSup, nrow, ncol, altHeadInf.shape[0])
                ajoutSimu = (se - ss) / duree
                # Watertable altitude value for considered day being interpolated
                s = ss + (ajoutSimu * pas)

                # deapth : topographical altitude  - watertable altitude
                #d = topoSimu[nrow][ncol] - s


                # Watertable altitude value for simulated period for reference simulation
                r = get_non_dry_cell_hds_value(
                    refHead, nrow, ncol, refHead.shape[0])

                WsRef = getWeightToSurface(topo_ref[nrow][ncol], r, dc, alpha)
                WsSimu = getWeightToSurface(topoSimu[nrow][ncol], s, dc, alpha)

                mWs = max(WsRef, WsSimu)
                sherrorsup += (mWs * (r-s)**2)
                smWs += mWs

    hErrorGlobal = math.sqrt(sherrorsup / smWs)

    with open(repo_simu + "/" + simu_name + '_Ref_' + ref_name + '_errorsresult_H_light.csv', 'w') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(['H Error'])
        writer.writerow([hErrorGlobal])

    return nbrowtot, nbcoltot



def compute_h_error_by_interpolation_cut(site_number, chronicle, approx, rate, ref, folder, time_step=1):
    # Get data for reference simulation
    # Path to repository
    ref_name = get_model_name(site_number, chronicle, approx, rate, ref=True)
    if folder is None:
        repo_ref = get_path_to_simulation_directory(site_number, chronicle, approx, rate, ref=True)
    else :
        site_name = get_site_name_from_site_number(site_number)
        repo_ref = folder + site_name + ref_name

    # contenu = os.listdir("/modflow/outputs/")
    # print("contenu", contenu)

    topo_ref = np.load(repo_ref + "/soil_surface_topo_"+ ref_name + ".npy")        
    # Watertable altitude
    #ref_hds = fpu.HeadFile(repo_ref + '/' + ref_name + '.hds')    
    
    
    if ref:
        # Get data for alternate simulation
        # Path to repository
        repo_simu = repo_ref 
        simu_name = ref_name
        # Topographical altitude
        topoSimu = topo_ref
        # Get heads values for simulation
        #simu_hds = ref_hds
        
    else:
        # Get data for alternate simulation
        # Path to repository
        simu_name = get_model_name(site_number, chronicle, approx, rate, ref=False)
        if folder is None:
            repo_simu = get_path_to_simulation_directory(site_number, chronicle, approx, rate, ref=False)
            
        else:
            repo_simu = folder + site_name + simu_name
        # Topographical altitude
        topoSimu = get_soil_surface_values_for_a_simulation(repo_simu, simu_name)
        # Get heads values for simulation

        ## test debug
        

        simu_hds = fpu.HeadFile(repo_simu + '/' + simu_name + '.hds')

    # Duration of simulated periods
    simu_times = simu_hds.get_times()

    # Initialisation of error computing variables
    smWs = 0
    sherrorsup = 0

    # For every day, we have to compare the reference and alternate simulations
    for day in range(startTime, endTime+1):
        # print(day)

        # Retrieve the watertable altitude values for the reference simulation for the considered day
        refHead = np.load(repo_ref + "/HDS/hds_" + str(day) + ".npy")
        #ref_hds.get_data(kstpkper=(0, day))

        # Compute the number of the simulated period of the alternate simulation for the considered day
        nbPeriod = 0
        while (simu_times[nbPeriod] < day+1) and (nbPeriod < len(simu_times)):
            nbPeriod += 1
        # print("nbPeriod : " + str(nbPeriod))

        # When the considered day match with the duration of the corresponding simulated period
        # We retrieve the value of the watertable altitude from the corresponding matrix
        if math.isclose(simu_times[nbPeriod], day+1, rel_tol=1e-3):
            altHeadSup = simu_hds.get_data(kstpkper=(time_step-1, nbPeriod))
            altHeadInf = altHeadSup
            duree = int(simu_times[nbPeriod])
            pas = 0
        # Otherwise, we have to interpolate the watertable altitude value for the considered day
        else:
            # The considered day is situated between the simulated period number 'nbPeriod-1' and number 'nbPeriod'
            altHeadSup = simu_hds.get_data(kstpkper=(time_step-1, nbPeriod))
            altHeadInf = simu_hds.get_data(kstpkper=(time_step-1, nbPeriod-1))
            duree = int(simu_times[nbPeriod] - simu_times[nbPeriod-1])
            pas = day - simu_times[nbPeriod-1]

        nbrowtot = altHeadInf.shape[1]
        nbcoltot = altHeadInf.shape[2]

        # We want to store the presence of cells part of a flood episode
        # We go through all the cells of the matrix representing the study site
        for nrow in range(nbrowtot):
            for ncol in range(nbcoltot):
                # Watertable altitude value for simulated period with duration lower than considered day
                ss = get_non_dry_cell_hds_value(
                    altHeadInf, nrow, ncol, altHeadInf.shape[0])
                # Watertable altitude value for simulated period with duration higher than considered day
                se = get_non_dry_cell_hds_value(
                    altHeadSup, nrow, ncol, altHeadInf.shape[0])
                ajoutSimu = (se - ss) / duree
                # Watertable altitude value for considered day being interpolated
                s = ss + (ajoutSimu * pas)

                # deapth : topographical altitude  - watertable altitude
                #d = topoSimu[nrow][ncol] - s


                # Watertable altitude value for simulated period for reference simulation
                r = get_non_dry_cell_hds_value(
                    refHead, nrow, ncol, refHead.shape[0])

                WsRef = getWeightToSurface(topo_ref[nrow][ncol], r, dc, alpha)
                WsSimu = getWeightToSurface(topoSimu[nrow][ncol], s, dc, alpha)

                mWs = max(WsRef, WsSimu)
                sherrorsup += (mWs * (r-s)**2)
                smWs += mWs

    hErrorGlobal = math.sqrt(sherrorsup / smWs)

    with open(repo_simu + "/" + simu_name + '_Ref_' + ref_name + '_errorsresult_H_light_cut.csv', 'w') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(['H Error'])
        writer.writerow([hErrorGlobal])

    return nbrowtot, nbcoltot



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-chr", "--chronicle", type=int, required=True)
    parser.add_argument("-approx", "--approximation", type=int, required=False)
    parser.add_argument("-rate", "--rate", type=float, required=False)
    parser.add_argument("-site", "--sitenumber", type=int, required=True)
    parser.add_argument("-ref", "--ref", action='store_true')
    # For evolutivity
    parser.add_argument("-step", "--step", type=int, required=False)
    parser.add_argument("-f", "--folder", type=str, required=False)
    parser.add_argument("-topo", "--topo", action='store_true')


    args = parser.parse_args()

    approx = args.approximation
    chronicle = args.chronicle
    site_number = args.sitenumber
    rate = args.rate
    ref = args.ref
    folder = args.folder
    topo = args.topo
    # dc = args.criticaldepth
    # alpha = args.alpha

    if topo:
        topo_file(site_number, chronicle, approx, rate, ref=True, folder=folder)
    else:
        nbrowtot, nbcoltot = compute_h_error_by_interpolation(site_number, chronicle, approx, rate, ref=False, folder=folder)
        print(nbrowtot, nbcoltot)
    #cut_hds_ref_file_into_one_day_hds_files(site_number, chronicle, approx, rate, folder)
    # topo_file(site_number, chronicle, approx, rate, ref, folder)