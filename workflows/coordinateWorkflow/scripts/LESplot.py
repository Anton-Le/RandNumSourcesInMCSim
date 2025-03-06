import numpy as np
from matplotlib import pyplot as plt
import os, shutil, sys, re

from scipy import stats

from scipy.spatial import cKDTree
from scipy import spatial

# function definitions
def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    deltas = nodes - node
    dist = np.linalg.norm(deltas, axis=1)
    min_idx = np.argmin(dist)
    return nodes[min_idx], dist[min_idx], deltas[min_idx][1]/deltas[min_idx][0]
def pointWithinBounds(point, maxBound, minBound):
    return ( point[0] >= minBound[0] and point[0] <= maxBound[0]) and ( point[1] >= minBound[1] and point[1] <= maxBound[1])
def circleWithinBounds(point, radius, maxBound, minBound):
    isWithinXboundaries = ( point[0]-radius >= minBound[0] and point[0]+radius <= maxBound[0])
    isWithinYboundaries = ( point[1]-radius >= minBound[1] and point[1]+radius <= maxBound[1])
    return isWithinXboundaries and isWithinYboundaries
def largestCircleRadius(points):
    """
    Returns the radius of the largest circle whose
    centre is still within the domain and 
    """
    if points.shape[0] >= 4:
        vor = spatial.Voronoi(points)
        max_d = 0
        max_v = None
        for v in vor.vertices:
            _, d, _ = closest_node(v, points)
            if d > max_d and circleWithinBounds(v, d, [1,1],[0,0]): #pointWithinBounds(v, [1,1], [0,0]):
                max_d = d
                max_v = v
    return max_d
def LESStatistics(dataset:np.array, numGroups:int):
    """
    This function splits the dataset of points into numGroups of
    equally sized groups and computes the discrepancy for each of the groups.
    """
    partitionedData = np.split( dataset, numGroups)
    discrepancies = np.zeros(numGroups, dtype=float)
    for subset, idx in zip(partitionedData, range(numGroups)):
        discrepancies[idx] = largestCircleRadius(subset)
    return discrepancies


figsize=(10,10)

#script parameters
baseDataDir=os.path.abspath(os.path.curdir)
fileName=snakemake.input[0]
groups=snakemake.params.groups
imageFileName=snakemake.output[0]

#derived parameters
FQFN=baseDataDir+'/'+fileName

#data processing
data = np.loadtxt(FQFN, usecols=[0,1], delimiter=',')

qrng_double_les = LESStatistics( data, groups)

LESList = [qrng_double_les]
datasetLabels=["single"]

# visualisation

numDatasets = len(LESList)

maxLES = snakemake.params.maxY # max( [max(item) for item in LESList] )
minLES = snakemake.params.minY # min( [min(item) for item in LESList] )

fig = plt.figure(figsize=(10, numDatasets *10))
gs = fig.add_gridspec(numDatasets, 2, width_ratios=(4,1), left=0.1, right=0.9, wspace=0.05, hspace=0.1)


scatterPlotAxes = [ fig.add_subplot(gs[idx,0]) for idx in range(numDatasets) ]
histogramAxes = [ fig.add_subplot(gs[idx,1], sharey=scatterPlotAxes[idx]) for idx in range(numDatasets) ]

for (ax, ax_histy, idx ) in zip(scatterPlotAxes, histogramAxes, range(numDatasets) ):
    ax.plot(np.arange(len(LESList[idx])), LESList[idx], 'bo', ms=4)
    ax_histy.hist(LESList[idx], bins=int( np.ceil(1 + np.log2( LESList[idx].shape[0] )) ), alpha=0.2, orientation='horizontal')
    ax_histy.tick_params(axis="y", labelleft=False)
    
    ax_histy.set_xlabel("Event count")
    ax.set_ylim([minLES, maxLES])
    ax.set_ylabel("Discrepancy", fontsize=14)
    ax.set_xlabel("Batch id.", fontsize=14)
    ax.set_title(datasetLabels[idx])

fig.savefig(imageFileName)

# print statistics information
print("Largest Empty Sphere statistics")
for (idx, dataset) in zip(range(len(LESList)), LESList):
    print(datasetLabels[idx])
    print("Mean: ", np.average(dataset), " variance: ", np.var(dataset, ddof=1))

#save computed discrepancies and their statistics
np.savetxt(snakemake.output[1], qrng_double_les)
np.savetxt(snakemake.output[2], np.array([np.average(qrng_double_les), np.var(qrng_double_les, ddof=1)]) )