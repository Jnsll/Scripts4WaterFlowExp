## Execution time

import math
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from numpy import median


def createExecTimePlotOfIndicatorFromCSVGlobalFile(indicator, rech=False, chronicle=False, sitename="Agon-Coutainville", refValues=True, log=False):
    repo = "/DATA/These/Projects/Model/app/exps/"
    name = "Exps_" + indicator + "_Indicator_" + sitename
    if refValues:
        name+= "_withref"
    if chronicle :
        name+= "_chronicle"
    if rech :
        name += "_rech"

    file = repo + name + ".csv"
    dfp = pd.read_csv(repo + file, sep=",")
    
    suffix = ""
    
    if rech:
        number_colors = 11
        dfp = dfp.reindex([10,0,1,2,3,4,5,6,7,8,9])
        approximation = "recharge threshold"
        suffix += "_rech"
    else:
        number_colors = 9
        dfp = dfp.reindex([8,6,1,0,3,7,4,5,2])
        approximation = "period duration"
    colors = sns.color_palette("hls", number_colors)


    if chronicle:
        suffix += "_chronicle"
    if refValues:
        suffix += "_withref"

    if log:
        for index, row in dfp.iterrows():
            print(row["H Error"], index)
        if index != 10:
            print(index)
            dfp["H Error"][index] = math.log(row["H Error"])
        else:
            dfp["H Error"][index] = -float('Inf')
            dfp["H Error"][index] = math.log(dfp["H Error"][index])

        suffix+= "_log"   


    a = sns.relplot(x="Execution Time", y=indicator + " Error", hue="Approximation", palette=colors[::-1],data=dfp)
    a.set(xlabel='Execution Time (s)', ylabel=indicator + ' Indicator (m)')
    plt.plot(dfp['Execution Time'], dfp[indicator + " Error"], linewidth=2, alpha=0.2)
    a.fig.suptitle("Evolution of " + indicator + " indicator value according to the execution time \n Approximation with " + approximation + "\n" + sitename + "site")
    plt.subplots_adjust(top=0.8)
    
    max_exec_time = dfp[dfp[indicator + " Error"]!=0]["Execution Time"]
    H_threshold = 0.1
    plt.plot([0, max_exec_time], [H_threshold, H_threshold], linewidth=2, alpha=0.7, color="Red", dashes=[6, 2])  
    
    a.savefig(repo + 'Plot_' + indicator + '_Indicator_ExecTime_' + sitename + suffix + '.png')

def createApproxPlotOfIndicatorFromCSVGlobalFile(indicator, rech=False, chronicle=False, sitename="Agon-Coutainville", refValues=True, log=False):
    repo = "/DATA/These/Projects/Model/app/exps/"
    name = "Exps_" + indicator + "_Indicator_" + sitename
    if refValues:
        name+= "_withref"
    if chronicle :
        name+= "_chronicle"
    if rech :
        name += "_rech"

    file = repo + name + ".csv"
    dfp = pd.read_csv(repo + file, sep=",")
    
    suffix = ""
    
    if rech:
        number_colors = 11
        dfp = dfp.reindex([10,0,1,2,3,4,5,6,7,8,9])
        approximation = "recharge threshold"
        suffix += "_rech"
    else:
        number_colors = 9
        dfp = dfp.reindex([8,6,1,0,3,7,4,5,2])
        approximation = "period duration"
    colors = sns.color_palette("hls", number_colors)


    if chronicle:
        suffix += "_chronicle"
    if refValues:
        suffix += "_withref"

    if log:
        for index, row in dfp.iterrows():
            print(row["H Error"], index)
        if index != 10:
            print(index)
            dfp["H Error"][index] = math.log(row["H Error"])
        else:
            dfp["H Error"][index] = -float('Inf')
            dfp["H Error"][index] = math.log(dfp["H Error"][index])

        suffix+= "_log"   


    a = sns.relplot(x="Approximation", y=indicator + " Error", hue="Approximation", palette=colors[::-1],data=dfp)
    a.set(xlabel=approximation, ylabel=indicator + ' Indicator (m)')
    plt.plot(dfp['Approximation'], dfp[indicator + " Error"], linewidth=2, alpha=0.2)
    a.fig.suptitle("Evolution of " + indicator + " indicator value according to the approximation rate \n Approximation with " + approximation + "\n" + sitename + "site")
    plt.subplots_adjust(top=0.8)
    
    max_exec_time = dfp[dfp[indicator + " Error"]!=0]["Execution Time"]
    H_threshold = 0.1
    plt.plot([0, max_exec_time], [H_threshold, H_threshold], linewidth=2, alpha=0.7, color="Red", dashes=[6, 2])  
    
    a.savefig(repo + 'Plot_' + indicator + '_Indicator_ApproxRate_' + sitename + suffix + '.png')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-ind", "--indicator", type=int, required=True)
    parser.add_argument("-site", "--site", help= "2: Agon-Coutainville or 3:Saint-Germain-Sur-Ay", type=str, required=True)
    parser.add_argument("-ref", "--reference", action='store_true')
    parser.add_argument("-chr", "--chronicle", action='store_true')
    parser.add_argument("-rech", "--recharge", action='store_true')
    parser.add_argument("-log", "--log", action='store_true')
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
    log = args.log
    createExecTimePlotOfIndicatorFromCSVGlobalFile(indicator, rech=rech, chronicle=chronicle, sitename=sitename, refValues=refValues, log=log)
