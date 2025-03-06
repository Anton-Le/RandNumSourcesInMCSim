import numpy as np
from matplotlib import pyplot as plt
import os, shutil, sys, re

from scipy import stats

from scipy.spatial import cKDTree

# function definitions
def discrepancyStatistics(dataset:np.array, numGroups:int):
    """
    This function splits the dataset of points into numGroups of
    equally sized groups and computes the discrepancy for each of the groups.
    """
    partitionedData = np.split( dataset, numGroups)
    discrepancies = np.zeros(numGroups, dtype=float)
    
    PE_Q = np.zeros(numGroups, dtype=float)
    avg_discrepancies = np.zeros(numGroups, dtype=float)
    CoV = np.zeros(numGroups, dtype=float)
    Rm = np.zeros(numGroups, dtype=float)

    for subset, idx in zip(partitionedData, range(numGroups)):
        #discrepancies[idx] = discrepancy(subset)
        
        # Create a KDTree from the points
        kdtree = cKDTree(subset)

        # Find the nearest neighbor for each point in the set
        nearest_neighbors = kdtree.query(subset, k=2)  # k=2 to exclude the point itself

        # Initialize a variable to store the largest distance
        largest_distance = 0.0

        smallest_distance = 10.0

        avg_distance = 0.0
        CV = 0.0

        Potential_energy_quality = 0.0
        n = len(subset)
        theta=3/n
    
        # Query for nearest neighbors in the entire dataset, excluding the point itself
        for point in subset:
            distances, indices = kdtree.query(point, k=2)  # k=2 to exclude the point itself
            max_distance = distances[1]  # The distance to the second nearest neighbor
            largest_distance = max(largest_distance, max_distance)

            smallest_distance = min(smallest_distance, max_distance)

            avg_distance = avg_distance + max_distance
            Potential_energy_quality = Potential_energy_quality + (1 - theta/(theta + max_distance**2))

        avg_distance = avg_distance/n

        discrepancies[idx] = largest_distance

        PE_Q[idx] = Potential_energy_quality/n
        Rm[idx] = largest_distance/smallest_distance

        for point in subset:
            distances, indices = kdtree.query(point, k=2)  # k=2 to exclude the point itself
            max_distance = distances[1]  # The distance to the second nearest neighbor
            CV = CV + (max_distance - avg_distance)**2
        CV = np.sqrt(1/n*CV)/avg_distance

        avg_discrepancies[idx] = avg_distance
        CoV[idx] = CV

    return discrepancies, PE_Q, Rm, avg_discrepancies, CoV

figsize=(10,10)

#script parameters
baseDataDir=os.path.abspath(os.path.curdir)
fileName=snakemake.input[0]
groups=snakemake.params.groups
imageFileName=snakemake.output[0]

# print info
print(baseDataDir)
print(fileName)
print(groups)
print(imageFileName)
#derived parameters
FQFN=baseDataDir+'/'+fileName

print(FQFN)

#data processing
data = np.loadtxt(FQFN, usecols=[0,1], delimiter=',')
#discrepancies = discrepancyStatistics( data, groups)
discrepancies, PE_Q, Rm, avg_discrepancies, CoV = discrepancyStatistics( data, groups)

discrepancyList = [discrepancies]
datasetLabels=[snakemake.wildcards.precision]

# visualisation

numDatasets = len(discrepancyList)

maxDiscrepancy = snakemake.params.maxY # max( [max(item) for item in discrepancyList] )
minDiscrepancy = snakemake.params.minY # min( [min(item) for item in discrepancyList] )

fig = plt.figure(figsize=(10, numDatasets *10))
gs = fig.add_gridspec(numDatasets, 2, width_ratios=(4,1), left=0.1, right=0.9, wspace=0.05, hspace=0.1)


scatterPlotAxes = [ fig.add_subplot(gs[idx,0]) for idx in range(numDatasets) ]
histogramAxes = [ fig.add_subplot(gs[idx,1], sharey=scatterPlotAxes[idx]) for idx in range(numDatasets) ]

for (ax, ax_histy, idx ) in zip(scatterPlotAxes, histogramAxes, range(numDatasets) ):
    ax.plot(np.arange(len(discrepancyList[idx])), discrepancyList[idx], 'bo', ms=4)
    ax_histy.hist(discrepancyList[idx], bins=int( np.ceil(1 + np.log2( discrepancyList[idx].shape[0] )) ), alpha=0.2, orientation='horizontal')
    ax_histy.tick_params(axis="y", labelleft=False)
    
    ax_histy.set_xlabel("Event count")
    ax.set_ylim([minDiscrepancy, maxDiscrepancy])
    ax.set_ylabel("Discrepancy", fontsize=14)
    ax.set_xlabel("Batch id.", fontsize=14)
    ax.set_title(datasetLabels[idx])

fig.savefig(imageFileName)

# print statistics information
print("Nearest neighbour distance statistics")
for (idx, dataset) in zip(range(len(discrepancyList)), discrepancyList):
    print(datasetLabels[idx])
    print("Mean: ", np.average(dataset), " variance: ", np.var(dataset, ddof=1))

#save computed discrepancies and their statistics
np.savetxt(snakemake.output[1], discrepancies)
np.savetxt(snakemake.output[2], np.array([np.average(discrepancies), np.var(discrepancies, ddof=1)]) )
np.savetxt(snakemake.output[3], PE_Q)
np.savetxt(snakemake.output[4], Rm)
np.savetxt(snakemake.output[5], avg_discrepancies)
np.savetxt(snakemake.output[6], CoV)
