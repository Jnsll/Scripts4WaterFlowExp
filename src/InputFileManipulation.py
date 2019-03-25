import os
import pandas as pd
import argparse

def extractDfFromRefInputFile(inputName):
    inputRepo = "resources"
    inputPath = os.path.join(
        '/'.join(os.getcwd().split('/')[:-1]), inputRepo, inputName)
    df = pd.read_csv(inputPath, sep='\t')
    return df

def getNumberOfLinesToReduceInLoop(prd):
    if (prd == "semester"):
        period = [183, 182]
    elif (prd == "monthly"):
        period = [30]
    elif (prd == "weekly"):
        period = [7]
    elif (prd == "daily"):
        period = [1]
    elif (prd == "init"):
        period = [1]
    return period

# nbTotalLines = len(df.index)


def getIndexesToRemain(nbTotalLines, period):
    ligne = 0
    indexes = [0]
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


def modifyValuesInLinesToRemain(df, indexes, nbTotal):
    for i in range(len(indexes)):
        if i < (len(indexes)-1):
            for z in range(indexes[i], indexes[i+1]):
                df.iat[indexes[i], 4] += df['rech'][z]
                #df.iat[indexes[i], 5] += df['chd_sea'][z]
            #df.iat[indexes[i], 5] = df.iat[indexes[i], 5] / \
             #   (indexes[i+1]-indexes[i])
            df.iat[indexes[i], 4] = df.iat[indexes[i], 4] / \
                (indexes[i+1]-indexes[i])
            df.iat[indexes[i], 1] = indexes[i+1]-indexes[i]
        else:
            for z in range(indexes[i], nbTotal):
                df.iat[indexes[i], 4] += df['rech'][z]
                #df.iat[indexes[i], 5] += df['chd_sea'][z]
            #df.iat[indexes[i], 5] = df.iat[indexes[i], 5]/(nbTotal-indexes[i])
            df.iat[indexes[i], 4] = df.iat[indexes[i], 4]/(nbTotal-indexes[i])
            df.iat[indexes[i], 1] = nbTotal-indexes[i]
    return df

def modifyValuesInLinesToRemainWithSeaLvl(df, indexes, nbTotal):
    for i in range(len(indexes)):
        if i < (len(indexes)-1):
            for z in range(indexes[i], indexes[i+1]):
                df.iat[indexes[i], 4] += df['rech'][z]
                df.iat[indexes[i], 5] += df['chd_sea'][z]
            df.iat[indexes[i], 5] = df.iat[indexes[i], 5] / \
                (indexes[i+1]-indexes[i])
            df.iat[indexes[i], 4] = df.iat[indexes[i], 4] / \
                (indexes[i+1]-indexes[i])
            df.iat[indexes[i], 1] = indexes[i+1]-indexes[i]
        else:
            for z in range(indexes[i], nbTotal):
                df.iat[indexes[i], 4] += df['rech'][z]
                df.iat[indexes[i], 5] += df['chd_sea'][z]
            df.iat[indexes[i], 5] = df.iat[indexes[i], 5]/(nbTotal-indexes[i])
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
    df["time_step"] = timestepValue
    return df


def changePeriodValue(df, periodValue):
    df["sp_length"] = periodValue
    return df


def addVariationToSeaLevel(df, seaLvlVariationValue):
    df["chd_sea"] += seaLvlVariationValue
    return df

def writeInputFile(modelname, df):
    outputname = "input_file_" + modelname + ".txt"
    filepath = os.path.join('/'.join(os.getcwd().split('/')[:-1]), "resources", outputname)
    df.to_csv(filepath, sep="\t", index=False)
    print(filepath)
    return filepath

def manipulateInputFile(inputName, prd, timestepValue, periodValue, seaValue):
    df = extractDfFromRefInputFile(inputName)
    if (prd != "init") and (prd is not None):
        period = getNumberOfLinesToReduceInLoop(prd)
        indexes = getIndexesToRemain(len(df.index), period)
        df = modifyValuesInLinesToRemain(df, indexes, len(df.index))
        df = removeLinesExceptThoseWithFollowingIndexes(df, indexes)
    if timestepValue is not None:
        df = changeTimestepValue(df, timestepValue)
    if periodValue is not None:
        df = changePeriodValue(df, periodValue)
    if seaValue is not None:
        df = addVariationToSeaLevel(df, seaValue)
    return df

def writeInputFileAfterManipulation(modelname, inputName, prd, timestepValue, periodValue, seaValue):
    df = manipulateInputFile(inputName, prd, timestepValue, periodValue, seaValue)
    filepath = writeInputFile(modelname, df)
    return filepath

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
    parser.add_argument(
        "-p", "--period", help="daily, weekly, monthly, semester or init", type=str, required=False)
    parser.add_argument("-ts", "--timestep", type=int, required=False)
    parser.add_argument("-sp", "--splength", type=int, required=False)
    parser.add_argument("-sea", "--sealvlchange", type=int, required=False)
    parser.add_argument("-m", "--modelname", type=str, required=True)
    parser.add_argument(
        "-i", "--inputfile", help="Name of the input file to take as reference", type=str, required=False)
    args = parser.parse_args()

    # Parameters
    if args.inputfile is None:
        inputName = "input_file.txt"
    else:
        inputName = args.inputfile
    timestepValue = args.timestep
    periodValue = args.splength
    modelname = args.modelname
    prd = args.period
    seaValue = args.sealvlchange

    writeInputFileAfterManipulation(modelname, inputName, prd, timestepValue, periodValue, seaValue)
