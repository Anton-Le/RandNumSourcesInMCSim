# README - coordinate workflow

## Workflow files

1. `Snakefile_QRNG.smk`
2. `Snakefile_PRNG.smk`
3. `coord_QRNG_snakemake.yaml`
4. `coord_PRNG_snakemake.yaml`

## Description

This workflow encompasses compilation of the code, generation of points on [0,1]^d
and aggreagation as well as visualisation of the statistics of the maximal nearest-neighbour distance
and the largest empty sphere radius.

**Caveat**: Currently the restirction of the evaluation to `d=2` applies.

The QRNG workflow performs, in general, the following steps:
1. Copy the source code directory tree from the directory provided as `repoPath` into
a directory in `./src` indexed by the `entropy` index (file number/seed number).
2. Create a directory in `./build` indexed by the `entropy` index and the chosen `precision` (**single, double only**),
substitute the file name provided as a source of entropy into the proxy functions in the code,
run within it `cmake` with appropriate compiler flags and compile the code using `make`.
3. Copy the resulting binary into the `./bin` directory, appending to the name the `entropy` index and `precision`
used.
4. Execute the code for the required number of samples and dimension in a directory in `./run`, to generate a `coordinates.lst` file
containing the coordinates of each sample point in [0,1]^d.
5. Collect nearest-neighbour distances and/or largest empty sphere radii of the generated data and for each `entropy` source
visualise the distribution of these in `./vis`.
6. Visualise - via a boxplot - the statistics of the means and variances of the LES radii.

QRNG and PRNG workflows differ, in whether a substitution is performed in step 3 (QRNG) or whether the `initSeed` parameter is passed
to the executable (PRNG) and whether the latter is compiled for QRNG file-based usage or PRNG usage.

*Note that usage of a QRNG device can be achieved by setting* `-DQRNG_USE_DEVICE=1` *in the workflow file, in which case the*
`entropy` *entries serve only as enumerators of cases to run and the values (seeds or file names) are de-facto ignored*.


## Running a workflow

### Prerequisites
Some rules of the workflow execute Python [scripts](.scripts) and hence require that the packages used
in those scripts are present in the Python environment in which Snakemake is installed.

These are as follows:
1. NumPy
2. SciPy (specifically the `stats` and `spatial` packages)
3. Matplotlib


To compile the C++ code the requirements for it have to be fulfilled, too:
1. Cmake (>3.20)
2. PCAP (optional)
3. QRNG library (optional).
4. Boost.

### Execution

1. To run _a_ workflow first create a directory within which the workflow will be executed.
2. Copy the appropriate snakefile, YAML configuration file and the entire `scripts` directory into the
newly created directory.
3. Rename the copied Snakefile, e.g. `Snakefile_PRNG.smk`, into `Snakefile`.
4. Load the Python virtual environment in which Snakemake has been installed.
5. Execute `snakemake --cores=1 -s Snakefile.smk vis/lesMeanValueStatistics.png` or `snakemake --cores=1 -s Snakefile.smk  vis/nnMeanValueStatistics.png`,
to collect statistics for 10 repetitions (each with a different entropy file/rng seed).

The generation of one (LES or NN) plot in step 5 will run through the entire process of compilation and point generation, but will only
execute the data aggregation and visualisation step (step 5 in the [##Description]) for the respective statistic (LES or NN).
Re-running Snakemake with the other plot afterwards will only run the rules necessary to aggregate and visualise the data.

To use more than 10 different RNG seeds or binary entropy files modify the appropriate configuration file and expand the list of seeds and 
change the upper limit of the range in line 4 of the associated Snakefile to cover the new range.

## Remarks

In the current state the PRNG and QRNG workflows are essentially duplicates and will have to be modularised and generalised
in the future, to avoid copy-paste errors.
