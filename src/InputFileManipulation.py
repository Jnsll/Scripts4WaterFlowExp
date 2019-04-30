import os
import pandas as pd
import argparse
import numpy as np

def extractDfFromRefInputFile(inputName):
    inputRepo = "data"
    inputPath = os.path.join('/'.join(os.getcwd().split('/')), inputRepo, inputName)
    df = pd.read_csv(inputPath, sep='\t', dtype=np.float64)
    return df

def getNumberOfLinesToReduceInLoop(prd):
    if (prd == "semester"):
        period = [182, 183, 182, 183, 183, 183, 183, 182]
    elif (prd== "trimester"):
        period = [91, 91, 91, 92, 91, 91, 91, 92, 92, 91, 91, 92, 91, 91, 91, 92]
    elif (prd == "monthly"):
        period = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        period.extend([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
        period.extend([31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
        period.extend([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    elif (prd == "weekly"):
        period = [7] * 52
        period.append(1)
        period.extend([7] * 52)
        period.append(1)
        period.extend([7]*52)
        period.append(2)
        period.extend([7] * 52)
        period.append(1)
    elif (prd == "twodays"):
        period = [2]
    elif (prd == "daily"):
        period = [1]
    elif (prd == "init"):
        period = [1]
    elif (prd == "year"):
        period = [365, 365, 366, 365]
    elif (prd == "twoyear"):
        period = [730, 731]
    elif (prd == "tenyear"):
        period = [3652, 3653]
    elif (prd == "tendays"):
        period = [10]
    elif (prd == "hundreddays"):
        period = [100]
    return period

# nbTotalLines = len(df.index)

def getIndexesToRemain(nbTotalLines, period):
    ligne = 1
    indexes = [0,1]
    taille_period = len(period)
    i = 0
    while (ligne < nbTotalLines-1):
        ligne += period[i]
        if (ligne < nbTotalLines-1):
            indexes.append(ligne)
        if (i == (taille_period-1)):
            i = 0
        else:
            i += 1
    return indexes

def getIndexesToRemainForRechThreshold(df, nbTotalLines, rechThr):
    indexes = [0,1]
    rechCompt = 0
    for ligne in range(1, nbTotalLines):
        rl = df['rech'][ligne]
        predCompt = rechCompt + rl
        if (ligne != nbTotalLines-1):
            if (predCompt > rechThr):
                if (predCompt - rechThr) < (rechThr - rechCompt):
                    indexes.append(ligne+1) # the line number ligne is included in perforation. The following line is the line to remain
                else:
                    indexes.append(ligne) # the line number ligne is not included in perforation and is the next line to remain
                rechCompt = 0
            else:
                rechCompt = predCompt
        else: # Last line
            if (predCompt > rechThr):
                if (predCompt - rechThr) >= (rechThr - rechCompt):
                    indexes.append(ligne)
    return indexes


def modifyValuesInLinesToRemain(df, indexes, nbTotal):
    df.loc[0,'rech'] = float(df['rech'].mean()) #df['rech'][0]
    print("mean init for inputfile: " + str(float(df['rech'].mean())))
    #print(df.loc[0,'rech'])
    for i in range(1, len(indexes)): #Initialisation has to remain so strating at period number 1. Initialisation period is number 0.
        if i < (len(indexes)-1):
            for z in range(indexes[i]+1, indexes[i+1]): #+1 pour ne pas compter la valeur df.iat[indexes[i]] une deuxième fois
                df.iat[indexes[i], 4] += df['rech'][z]
            df.iat[indexes[i], 4] = df.iat[indexes[i], 4] / \
                (indexes[i+1]-indexes[i])
            df.iat[indexes[i], 1] = indexes[i+1]-indexes[i]
        else:
            for z in range(indexes[i]+1, nbTotal):
                df.iat[indexes[i], 4] += df['rech'][z]
            df.iat[indexes[i], 4] = df.iat[indexes[i], 4]/(nbTotal-indexes[i])
            df.iat[indexes[i], 1] = nbTotal-indexes[i]
    return df

def modifyValuesInLinesToRemainWithSeaLvl(df, indexes, nbTotal):
    for i in range(len(indexes)):
        if i < (len(indexes)-1):
            for z in range(indexes[i]+1, indexes[i+1]):
                df.iat[indexes[i], 4] += df['rech'][z]
            df.iat[indexes[i], 4] = df.iat[indexes[i], 4] / \
                (indexes[i+1]-indexes[i])
            df.iat[indexes[i], 1] = indexes[i+1]-indexes[i]
        else:
            for z in range(indexes[i]+1, nbTotal):
                df.iat[indexes[i], 4] += df['rech'][z]
            df.iat[indexes[i], 4] = df.iat[indexes[i], 4]/(nbTotal-indexes[i])
            df.iat[indexes[i], 1] = nbTotal-indexes[i]
    return df


def removeLinesExceptThoseWithFollowingIndexes(df, indexes):
    index_remove = []
    for y in range(0, len(df.index)):
        index_remove.append(y)

    for z in indexes[:]:
        if z in index_remove:
            index_remove.remove(z)

    df.drop(df.index[index_remove], inplace=True)
    return df


def changeTimestepValue(df, timestepValue):
    df["time_step"][1:] = timestepValue #[1:] when initialisation was conserved
    return df


def changePeriodValue(df, periodValue):
    df["sp_length"] = periodValue
    return df


def writeInputFile(modelname, df):
    outputname = "input_file_" + modelname + ".txt"
    filepath = os.path.join('/'.join(os.getcwd().split('/')), "data", outputname)
    df.to_csv(filepath, sep="\t", index=False)
    return outputname

def manipulateInputFile(inputName, prd, timestepValue, rechThreshold): #, periodValue
    df = extractDfFromRefInputFile(inputName)

    if rechThreshold is not None:
        indexes = getIndexesToRemainForRechThreshold(df, len(df.index), rechThreshold)

    if (prd != "init") and (prd is not None):
        period = getNumberOfLinesToReduceInLoop(prd)
        indexes = getIndexesToRemain(len(df.index), period)

    df = modifyValuesInLinesToRemain(df, indexes, len(df.index))
    df = removeLinesExceptThoseWithFollowingIndexes(df, indexes)
    
    if timestepValue is not None:
        df = changeTimestepValue(df, timestepValue)

    return df

def writeInputFileAfterManipulation(modelname, prd, timestepValue, rechThreshold, inputName):
    print("inp :" + inputName)
    df = manipulateInputFile(inputName, prd, timestepValue, rechThreshold) #, periodValue
    outputname = writeInputFile(modelname, df)
    return outputname

if __name__ == '__main__':
    # if len(sys.argv) == 1:
    #     servername = 'localhost'
    #     port = 1883
    # elif len(sys.argv) == 2:
    #     servername = sys.argv[1]
    #     port = 1883
    # elif len(sys.argv) == 3:
    #     servername = sys.argv[1]
    #     port = sys.argv[2]
    # else:
    #     print(__doc__)

    ### Parser ###
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--period", type=str, required=False)
    parser.add_argument("-ts", "--timestep", type=int, required=False)
    parser.add_argument("-r", "--rechthreshold", type=int, required=False)
    parser.add_argument("-m", "--modelname", type=str, required=True)
    parser.add_argument("-i", "--inputfile", help="Name of the input file to take as reference", type=str, required=False)
    args = parser.parse_args()

    # Parameters
    if args.inputfile is None:
        inputName = "input_file.txt"
    else:
        inputName = args.inputfile
    timestepValue = args.timestep
    modelname = args.modelname
    rechThreshold = args.rechthreshold
    prd = args.period

    writeInputFileAfterManipulation(modelname, prd, timestepValue, rechThreshold, inputName)
