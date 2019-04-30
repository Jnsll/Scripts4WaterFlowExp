import re
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def getConcatDfFromExecutionTimeAndErrorDf(dfDiff, l, taille, exec_t, duration):
    a_exec_time = np.array([exec_t]*taille)
    df_exectime = pd.DataFrame(data=a_exec_time, index=l.index, columns=['Execution Time'])
    df_dur = pd.DataFrame(data=[duration]*taille, index=l.index, columns=['Period Duration'])
    ens = pd.concat([l, df_exectime, df_dur], axis=1)
    dfDiff = pd.concat([dfDiff, ens])
    return dfDiff

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
    else:
        exec_time = int(m.group(1))*24*60*60 + int(m.group(2))*60*60 + int(m.group(3))*60 + float(m.group(4))
    return int(exec_time)


mainRepo = "/DATA/These/Projects/Model/app/Agon-Coutainville/"
ref = "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1"

duration = 3562
repo10yearsCsv = "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodtenyear_Step1/model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodtenyear_Step1_Ref_model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_errorsresult.csv"

dfDiff = pd.DataFrame(columns=['Period Number', 'Simulated Time', 'MAE', 'MRE', 'RMSE', 'Execution Time', 'Period Duration'])
df10y = pd.read_csv(mainRepo + repo10yearsCsv, sep= ";")

times = df10y["Simulated Time"]

file = mainRepo + "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodtenyear_Step1/model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodtenyear_Step1.list" 



exec_time = getExecutionTimeFromListFile(file)
dfDiff = getConcatDfFromExecutionTimeAndErrorDf(dfDiff, df10y, len(times), exec_time, duration)

simulations = ["model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodtwoyear_Step1", "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodyear_Step1", "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodsemester_Step1", "model_time_0_geo_0_thick_1_K_86.4_Sy_0.1_Periodtrimester_Step1"]
durations = [730, 365, 182, 91]
for nb in range(len(simulations)):
    dfSimu = pd.read_csv(mainRepo + simulations[nb] + "/" + simulations[nb] + "_Ref_" + ref + "_errorsresult.csv", sep=";")
    z=pd.DataFrame()
    for t in times:
        l = dfSimu.loc[dfSimu["Simulated Time"] == t]
        z = pd.concat([z,l])
    exec_t = getExecutionTimeFromListFile(mainRepo + simulations[nb] + "/" + simulations[nb] + ".list")
    dfDiff = getConcatDfFromExecutionTimeAndErrorDf(dfDiff, z, len(times), exec_t, durations[nb])


dfDiff.to_csv(os.path.join('/'.join(mainRepo.split('/')[:-2]), "exps/") + "exps_sameNbDiff.csv")


################## Plotting ################


exps = pd.read_csv("/DATA/These/Projects/Model/app/exps/exps_sameNbDiff.csv", sep=",")
fig, axs = plt.subplots(ncols=2)

plt.xticks(rotation=45)
colors = sns.color_palette("hls", 5)
a = sns.boxplot(x="Period Duration", y="RMSE", data=exps, ax=axs[0], palette = colors)
a.set(xlabel='Period Duration (days)', ylabel='Different of \ngroundwater table roof values(m)')


b = sns.boxplot(x="Execution Time", y="RMSE", data=exps, ax=axs[1], palette=colors[::-1])
b.set(xlabel='Execution Time (seconds)', ylabel='Different of \ngroundwater table roof values(m)')

for ax in fig.axes:
    plt.sca(ax)
    plt.xticks(rotation=45)
fig.tight_layout()

fig.savefig('/DATA/These/Projects/Model/app/exps/'+'Boxplot_diffwatertable_accordingto_durationOfperiod_executionTime_SameNumberOfDiffs.png')
