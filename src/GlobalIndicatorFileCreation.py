import re
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse


##### Generating CSV File with indicators for all simulations #####


def getExecutionTimeFromListFile(file):
    with open(file,'r') as f:
        lines = f.readlines()
        if lines:
            beforelast_line = lines[-2]
    beforelast_line=beforelast_line.rstrip()

    m = re.search(r'\sElapsed run time:\s+(?:(\d*)?(?:\sDays,\s*))?(?:(\d*)(?:\s*Hours,\s*))?(?:(\d*)(?:\s*Minutes,\s*))?(\d*[.]*\d*)\sSeconds', beforelast_line)
    if (m.group(1) is None) and (m.group(2) is None) and (m.group(3) is None):
        exec_time = float(m.group(4))
    elif (m.group(1) is None) and (m.group(2) is None):
        exec_time = int(m.group(3))*60 + float(m.group(4))      
    elif (m.group(1) is None):
        exec_time = int(m.group(2))*60*60 + int(m.group(3))*60 + float(m.group(4))
        print(m.group(2), m.group(3), m.group(4), exec_time)
    else:
        exec_time = int(m.group(1))*24*60*60 + int(m.group(2))*60*60 + int(m.group(3))*60 + float(m.group(4))
    return int(exec_time)

def createGlobalCSVFile(indicator, rech=False, chronicle=False, sitename=2, refValues=True):
    """
    sitename : #Agon-Coutainville #Saint-Germain-Sur-Ay
    
    """
    
    mainRepo = "/DATA/These/Projects/Model/app/" + sitename + "/"
    
    ref = "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Perioddaily_Step1"

    
    if rech:
        modelnames = ["model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_RechThreshold0.0002",
                 "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_RechThreshold0.05",
                 "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_RechThreshold0.1",
                 "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_RechThreshold0.2",
                 "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_RechThreshold0.25",
                 "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_RechThreshold0.4",
                 "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_RechThreshold0.5",
                 "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_RechThreshold0.8",
                 "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_RechThreshold1.0",
                 "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Step1_RechThreshold2.0"
             ]
        approx = ["0.0002",
             "0.05",
             "0.1",
             "0.2",
             "0.25",
             "0.4",
             "0.5",
             "0.8",
             "1.0",
             "2.0"
             ]
        if refValues:
            approx.append("0")
        
    else :
        modelnames = ["model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodmonthly_Step1",
               "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodweekly_Step1",
               "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodtenyear_Step1",
               "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodtrimester_Step1",
              "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodyear_Step1",
              "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodtwoyear_Step1",
              "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodtwodays_Step1",
              "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodsemester_Step1"          
           ]
        approx = ["30",
          "7",
          "3652",
          "91",
          "365",
          "730",
          "2",
          "182"
          ]
        if refValues:
            approx.append("1")
   
    
    if refValues:
        modelnames.append("model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Perioddaily_Step1")
        
    if chronicle:
        ref += "_chronicle"
           
    dfglob = pd.DataFrame()

    for ind in range(len(modelnames)):
        simuRepo = modelnames[ind]
        if refValues and (ind == len(modelnames)-1):
            filename = modelnames[ind] + "_Ref_" + "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1" + "_errorsresult_" + indicator +".csv"
        else:
        
            filename = modelnames[ind] + "_Ref_" + ref + "_errorsresult_" + indicator + "_light.csv" 
        
        if chronicle :
            simuRepo += "_chronicle"
            filename = modelnames[ind]  + "_chronicle" + "_Ref_" + ref + "_errorsresult_" + indicator + "_light.csv" 

        df = pd.read_csv(mainRepo + simuRepo + "/" + filename, sep=";")

    
        taille = len(df.index)
        df_dur = pd.DataFrame(data=[approx[ind]]*taille, index=df.index, columns=['Approximation'])
        exec_time = getExecutionTimeFromListFile(mainRepo + simuRepo + "/" + simuRepo + ".list")
        df_time = pd.DataFrame(data=[exec_time]*taille, index=df.index, columns=['Execution Time'])
        df = pd.concat([df, df_dur, df_time], axis=1)
        dfglob = pd.concat([dfglob,df])

    output = "Exps_" + indicator + "_Indicator_" + sitename
    if refValues:
        output+= "_withref"
    if chronicle :
        output+= "_chronicle"
    if rech :
        output += "_rech"
    
    dfglob.to_csv("/DATA/These/Projects/Model/app/exps/" + output + ".csv")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-ind", "--indicator", type=int, required=True)
    parser.add_argument("-site", "--site", help= "2: Agon-Coutainville or 3:Saint-Germain-Sur-Ay", type=str, required=True)
    parser.add_argument("-ref", "--reference", action='store_true')
    parser.add_argument("-chr", "--chronicle", action='store_true')
    parser.add_argument("-rech", "--recharge", action='store_true')
    args = parser.parse_args()
    
    if args.site == 2:
        sitename = "Agon-Coutainville"
    elif args.site == 3:
        sitename = "Saint-Germain-Sur-Ay"
    chronicle = args.chronicle
    refValues = args.reference
    outliers = True
    rech=args.recharge
    indicator = args.indicator
    createGlobalCSVFile(indicator, rech=rech, chronicle=chronicle, sitename=sitename, refValues=refValues)