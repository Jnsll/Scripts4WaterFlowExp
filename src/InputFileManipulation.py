import pandas as pd 


#Parameters
repo = "/DATA/These/Projects/WaterFlowSimulationModel/"
filename = "input_file_model_test_withoutFirst_supprlignes1sur2.txt"


df = pd.read_csv(repo+filename, sep='\t')
dfbis = df


indexes = [0,183,365,548]
nbToKeep = 4
nbTotal = len(dfbis.index)

def modifyValuesInLinesToRemain(indexes, nbTotal):
    for i in range(len(indexes)):
        if i <(len(indexes)-1):
            for z in range(indexes[i], indexes[i+1]):
                dfbis.iat[indexes[i],4] += dfbis['rech'][z]
                dfbis.iat[indexes[i],5] += dfbis['chd_sea'][z]
            dfbis.iat[indexes[i],5] = dfbis.iat[indexes[i],5]/(indexes[i+1]-indexes[i])
            dfbis.iat[indexes[i], 1] = indexes[i+1]-indexes[i]
        else:
            for z in range(indexes[i], nbTotal):
                dfbis.iat[indexes[i],4] += dfbis['rech'][z]
                dfbis.iat[indexes[i],5] += dfbis['chd_sea'][z]
            dfbis.iat[indexes[i],5] = dfbis.iat[indexes[i],5]/(nbTotal-indexes[i])
            dfbis.iat[indexes[i],1] = nbTotal-indexes[i]
    return dfbis
          
        
def removeLinesExceptThoseWithFollowingIndexes(dfbis, indexes):
    index_remove = []
    for y in range (0,nbTotal):
        index_remove.append(y)
    

    for z in indexes[:]:
        if z in index_remove:
            index_remove.remove(z)


    dfbis.drop(df.index[index_remove], inplace=True)
    return dfbis




if __name__ == '__main__':
    dfbis = modifyValuesInLinesToRemain(indexes, nbTotal)
    dfbiss = removeLinesExceptThoseWithFollowingIndexes(dfbis, indexes)
    outputname = "input_file_model_test_withoutFirst_semestre_comb.txt"
    dfbiss.to_csv(repo + outputname, sep="\t", index=False)