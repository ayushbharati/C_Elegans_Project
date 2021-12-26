import pandas as pd
from pandas._libs.missing import NA
import re
import requests
import time
import ast
import numpy as np
import matplotlib.pyplot as plt

def getCellTypesRaw(wbID):
    url = 'https://wormbase.org//rest/widget/gene/' + str(wbID) + '/overview'
    # page = requests.get(URL)
    page = ''
    while page == '':
        try:
            page = requests.get(url)
            break
        except:
            print("Connection refused by the server..")
            time.sleep(5)
            print("Was a nice sleep, now let me continue...")
            continue
    regexSearch = re.search("Is expressed in (.+?)[\.]", str(page.content))
    if regexSearch:
        print(regexSearch.group(1))
        return regexSearch.group(1)
    else:
        # print("NA -- regex error")
        return NA

def strToList(cellTypeStr):
    cellTypeStr = cellTypeStr.replace("several structures, including","")
    cellTypeStr = cellTypeStr.replace(";",",")
    cellTypeStr = cellTypeStr.replace("and",",")
    cellTypeStr = cellTypeStr.replace(" ","")
    cellTypes = cellTypeStr.split(',')

    for i in cellTypes:
        if i == "":
            cellTypes.remove(i)
    print(cellTypes)
    return cellTypes

def makeBarGraph(graphData: pd.DataFrame, fname = ""):
    N = graphData.shape[1] #number of column
 
    up = graphData.loc['up']
    down = graphData.loc['down']
    ind = np.arange(N)  
    width = 0.35 
    
    fig = plt.subplots(figsize =(10, 7))
    p1 = plt.bar(ind, up, width, bottom = down)
    p2 = plt.bar(ind, down, width)
    
    plt.ylabel('Number of genes')
    plt.title('Number of genes up/down regulated by cell type')
    plt.xticks(ind, graphData.columns)
    plt.yticks(np.arange(0, graphData.loc['total'].max() + 100, 100))
    plt.legend((p1[0], p2[0]), ('Up regulated', 'Down Regulated'))
    
    plt.savefig(fname + "barGraph.png")

inputFile = "FRKvsWT_significant_upvsdowngenes_22June"
# # # #
# reading and converting data

data = pd.read_csv("./data/" + inputFile + ".csv")   # read input csv
data['cellTypesRaw'] = data['WormBaseID'].apply(getCellTypesRaw)    # get expression info from wormbase
data = data[data['cellTypesRaw'].notna()]   #remove genes that didnt find cell types
data = data.reset_index()   # reindex df
data['cellTypesList'] = data['cellTypesRaw'].apply(strToList) # format cell type string to a list format

data.to_csv("./" + inputFile + "cellTypesFilteredLists.csv") # save processed data

# if reading from csv
# data = pd.read_csv("./" + inputFile + "cellTypesFilteredLists.csv")
# data['cellTypesList'] = data['cellTypesList'].apply(lambda x: ast.literal_eval(x)) # convert list string to list type

# # # #
# graphing data
geneCellMap = {}    # empty dictionary

# group row indexes by cell type into a dictionary
# this data structure will help with organizing the data before graphing
for index, row in data.iterrows():
    for cellType in row['cellTypesList']:
        if cellType in geneCellMap:
            geneCellMap[cellType].append(index)
        else:
            geneCellMap[cellType] = [index]

# prune keys that have few elements
threshold = 50 # number of genes a cell type needs to have to be considered

original_size = len(geneCellMap)
geneCellMap = {key:val for key, val in geneCellMap.items() if len(val) > threshold}
pruned_size = len(geneCellMap)
print("original: ", original_size, " pruned", pruned_size)

# sort dictionary based on number of genes in cell type
geneCellMapSort = sorted(geneCellMap.items(), key = lambda kv: len(kv[1]), reverse=True)

# init dataframe
graphDf = pd.DataFrame(columns = [cellType[0] for cellType in geneCellMapSort], index = ['up','down'])
graphDf.to_csv("./" + inputFile + "cellTypesUp_DownRegulated.csv")

for cellType in geneCellMapSort:
    up = 0
    down = 0
    for geneIndex in cellType[1]:
        if data.iloc[geneIndex]['logFC'] > 0:
            # up reg
            up = up + 1
        else:
            # down reg
            down = down + 1
    graphDf.loc['up', cellType[0]] = up
    graphDf.loc['down', cellType[0]] = down
    graphDf.loc['total', cellType[0]] = up + down
print(graphDf)

makeBarGraph(graphDf)