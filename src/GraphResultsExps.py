## Execution time

import math
import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from numpy import median

mainAppRepo = os.path.dirname(os.path.abspath(__file__)) + '/'

def get_site_name_from_site_number(site_number):
    sites = pd.read_csv(mainAppRepo + 'data/study_sites.txt',
                        sep=',', header=0, index_col=0) #\\s+
    site_name = sites.index._data[site_number]
    return site_name


def createExecTimePlotOfIndicatorFromCSVGlobalFile(indicator, folder, site, chronicle, approx, log=False):
    
    repo = folder
    site_name = get_site_name_from_site_number(site)
    name = "Exps_" + indicator + "_Indicator_" + site_name + "_Chronicle"+ str(chronicle) + "_Approx" + str(approx)


    file = repo + site_name + "/" + name + ".csv"
    dfp = pd.read_csv(file, sep=",")
    
    suffix = ""
    
    if approx == 1:
        number_colors = 11
        #dfp = dfp.reindex([10,0,1,2,3,4,5,6,7,8,9])
        approximation = "recharge threshold"
        #suffix += "_rech"
    else:
        number_colors = 9
        #dfp = dfp.reindex([8,6,1,0,3,7,4,5,2])
        approximation = "period duration"
    colors = sns.color_palette("hls", number_colors)


    if log:
        for index, row in dfp.iterrows():
            #print(row["H Error"], index)
            if index != 8:
                #print(index)
               #dfp.loc[index,"H Error"] = math.log(row["H Error"])
                dfp.loc[index,"Execution Time"] = math.log(row["Execution Time"])
            else:
                #dfp.loc[index,"H Error"] = -float('Inf')
                #dfp.loc[index,"H Error"] = math.log(-float('Inf'))
                dfp.loc[index,"Execution Time"] = math.log(row["Execution Time"])
                #dfp.loc[index,"Execution Time"] = math.log(-float('Inf'))

        suffix+= "_log"   


    a = sns.relplot(x="Execution Time", y=indicator + " Error", hue="Approximation", palette=colors[::-1],data=dfp)
    a.set(xlabel='Execution Time (s)', ylabel=indicator + ' Indicator (m)')
    plt.plot(dfp['Execution Time'], dfp[indicator + " Error"], linewidth=2, alpha=0.2)
    a.fig.suptitle("Evolution of " + indicator + " indicator value according to the execution time \n Approximation with " + approximation + "\n" + site_name + " site")
    plt.subplots_adjust(top=0.8)
    
    max_exec_time = dfp[dfp[indicator + " Error"] == 0]["Execution Time"] #-float('Inf')
    print(max_exec_time)
    if log:
        H_threshold = math.log(0.1)
    else:
        H_threshold = 0.1
    H_threshold = 0.1
    plt.plot([0, max_exec_time], [H_threshold, H_threshold], linewidth=3, alpha=0.7, color="Red", dashes=[6, 2])  
    
    a.savefig(repo + site_name + "/" + 'Plot_' + indicator + '_Indicator_ExecTime_' + site_name + "_Chronicle"+ str(chronicle) + "_Approx" + str(approx) + suffix + '.png')


def createTwoApproxExecTimePlotOfIndicatorFromCSVGlobalFile(indicator, rech=False, chronicle=False, sitename="Agon-Coutainville", refValues=True, log=False):
    repo = "/DATA/These/Projects/Model/app/exps/"
    name = "Exps_" + indicator + "_Indicator_" + sitename
    if refValues:
        name+= "_withref"
    if chronicle :
        name+= "_chronicle"
    # if rech :
    #     name += "_rech"

    file = repo + name + ".csv"
    file_rech = repo + name + "_rech.csv"
    dfp = pd.read_csv(file, sep=",")
    dfp_rech = pd.read_csv(file_rech, sep=",")
    
    suffix = ""
    

    number_colors_rech = 11
    dfp_rech = dfp_rech.reindex([10,0,1,2,3,4,5,6,7,8,9])
    approximation = "recharge threshold"

    number_colors = 9
    dfp = dfp.reindex([8,6,1,0,3,7,4,5,2])
    approximation = "period duration"
    

    df_type = pd.DataFrame(data=["Duration"]*len(dfp.index), index=dfp.index, columns=['Approx Type'])
    dfp = pd.concat([dfp, df_type], axis=1)

    df_type_rech = pd.DataFrame(data=["Recharge"]*len(dfp_rech.index), index=dfp_rech.index, columns=['Approx Type'])
    dfp_rech = pd.concat([dfp_rech, df_type_rech], axis=1)

    data = pd.concat([dfp, dfp_rech])


    if chronicle:
        suffix += "_chronicle"
    if refValues:
        suffix += "_withref"

    suffix += "_test_double_curves"
    

    a = sns.relplot(x="Execution Time", y=indicator + " Error", hue='Approx Type', data=data)
    a.set(xlabel='Execution Time (s)', ylabel=indicator + ' Indicator (m)')
    plt.plot(dfp['Execution Time'], dfp[indicator + " Error"], linewidth=2, alpha=0.2)
    plt.plot(dfp_rech['Execution Time'], dfp_rech[indicator + " Error"], linewidth=2, alpha=0.2)
    a.fig.suptitle("Evolution of " + indicator + " indicator value according to the execution time \n Approximation with " + approximation + "\n" + sitename + " site")
    plt.subplots_adjust(top=0.8)
    
    # max_exec_time = dfp[dfp[indicator + " Error"] == 0]["Execution Time"] #-float('Inf')

    # if log:
    #     H_threshold = math.log(0.1)
    # else:
    #     H_threshold = 0.1
    # H_threshold = 0.1
    # plt.plot([0, max_exec_time], [H_threshold, H_threshold], linewidth=3, alpha=0.7, color="Red", dashes=[6, 2])  
    
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
                dfp.loc[index,"H Error"] = math.log(row["H Error"])
                dfp.loc[index,"Execution Time"] = math.log(row["Execution Time"])
            else:
                dfp.loc[index,"H Error"] = -float('Inf')
                dfp.loc[index,"H Error"] = math.log(row["H Error"])
                dfp.loc[index,"Execution Time"] = -float('Inf')
                dfp.loc[index,"Execution Time"] = math.log(row["Execution Time"])

        suffix+= "_log"   


    a = sns.relplot(x="Approximation", y=indicator + " Error", hue="Approximation", palette=colors[::-1],data=dfp)
    a.set(xlabel=approximation, ylabel=indicator + ' Indicator (m)')
    plt.plot(dfp['Approximation'], dfp[indicator + " Error"], linewidth=2, alpha=0.2)
    a.fig.suptitle("Evolution of " + indicator + " indicator value according to the approximation rate \n Approximation with " + approximation + "\n" + sitename + "site")
    plt.subplots_adjust(top=0.8)
    
    max_exec_time = dfp[dfp[indicator + " Error"]!=0]["Execution Time"]
    if log:
        H_threshold = math.log(0.1)
    else:
        H_threshold = 0.1
    plt.plot([0, max_exec_time], [H_threshold, H_threshold], linewidth=2, alpha=0.7, color="Red", dashes=[6, 2])  
    
    a.savefig(repo + 'Plot_' + indicator + '_Indicator_ApproxRate_' + sitename + suffix + '.png')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-ind", "--indicator", type=str, required=True)
    parser.add_argument("-site", "--site", type=int, help= "2: Agon-Coutainville or 3:Saint-Germain-Sur-Ay", required=True)
    parser.add_argument("-chr", "--chronicle", type=int)
    parser.add_argument("-approx", "--approximation", type=int)
    parser.add_argument("-log", "--log", action='store_true')
    parser.add_argument("-f", "--folder", type=str, required=True)

    args = parser.parse_args()
    
    site = args.site 
    chronicle = args.chronicle
    outliers = True
    approx=args.approximation
    indicator = args.indicator
    log = args.log
    folder = args.folder
    #createTwoApproxExecTimePlotOfIndicatorFromCSVGlobalFile(indicator, rech=rech, chronicle=chronicle, sitename=sitename, refValues=refValues, log=log)
    createExecTimePlotOfIndicatorFromCSVGlobalFile(indicator, folder, site, chronicle, approx, log=False)