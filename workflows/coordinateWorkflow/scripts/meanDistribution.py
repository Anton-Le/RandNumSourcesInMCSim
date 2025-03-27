import numpy as np
from matplotlib import pyplot as plt
import os, shutil, sys, re

# function definitions
def collectValues(fileNameList:list):
    """
    The function collects
    """
    meanValues = np.zeros( len(fileNameList) )
    variances = np.zeros( len(fileNameList) )
    for idx, filename in zip(range(len(fileNameList)), fileNameList) :
        print(filename)
        data = np.loadtxt(filename)
        meanValues[idx] = data[0]
        variances[idx] = data[1]
    return (meanValues, variances)

figsize=(10,10)

#script parameters
baseDataDir=os.path.abspath(os.path.curdir)
fileNameList=snakemake.input
imageFileName=snakemake.output[0]

#data processing
means, var = collectValues(fileNameList)
datasetSize = means.shape[0]
print(datasetSize)
assert datasetSize % 2 == 0, "Not comparing two cases!" 

reshapedMeans = means.reshape(2, datasetSize//2 )
reshapedVar = var.reshape(2, datasetSize//2)
# visualisation

fig, (axm, axv) = plt.subplots(1, 2, sharex=True)

#draw means
axm.boxplot(reshapedMeans.transpose())

#draw variances
axv.boxplot(reshapedVar.transpose())

# common options
axm.set_xlabel("precision", fontsize=14)
axm.set_ylabel("Mean", fontsize=14)
axv.set_ylabel("Var", fontsize=14)

axm.set_xticks([1,2])
axm.set_xticklabels(["single", "double"], fontsize=14)



fig.savefig(imageFileName)