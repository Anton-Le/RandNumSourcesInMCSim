# README - pi workflow

## Workflow files

1. `Snakefile_QRNG.smk`
2. `Snakefile_PRNG.smk`
3. `configuration_PRNG.yaml`
4. `configuration_QRNG.yaml`

## Description

This workflow encompasses compilation of the code and generation of $\pi$-samples that should then be aggregated using
workflow processing notebooks within this repository.

The workflow performs, in general, the following steps:
1. Copy the source code directory tree from the the provided `repoPath` into
a directory in `./src` indexed by the `entropy` index (file/seed number).
2. Create a directory in `./build` indexed by the `entropy` index and the chosen `precision` (**single/double only**),
substitute the file name provided as a source of entropy into the proxy functions in the code (if the workflow is for QRNGs),
run within it `cmake` with appropriate compiler flags and compile the code using `make`.
3. Copy the resulting binary into the `./bin` directory's subdirectory indexed by the entropy file index and the `precision`.
4. Copy an appropriate assay bash script  template from the `scripts` directory to the `run/entropyFileIdx/precision` directory.
5. Instantiate the assay script by filling it with parameters provided in the configuration file as well as the paths of the binary to
use for the simulation.
6. Collect the data by running the assay script.

The workflow is considered complete when for each provided entropy file or PRNG seed and each requested precision the file
`run/<entropyFileIndex>/<precision>/completed` is present.

QRNG and PRNG workflows differ, in whether a substitution is performed in step 3 (QRNG) or whether the `initSeed` parameter is passed
to the executable (PRNG) and whether the latter is compiled for QRNG file-based usage or PRNG usage.

*Note that usage of a QRNG device can be achieved by setting* `-DQRNG_USE_DEVICE=1` *in the workflow file, in which case the*
`entropy` *entries serve only as enumerators of cases to run and the values ( seeds or file names) are de-facto ignored*.
This requires the modification of the `runCompilation` rule's `-DCMAKE_CXX_FLAGS` parameters as described in the top-level
README. By default the QRNG workflow will use `/dev/random` as a substitute for the QRNG.

**Also note** that the `repoPath` is copied in its entirety. Hence the directory the workflow is being run from _should not_ be
a sub-directory of the repository. It is additionally advised to remove any large datasets from the repository prior to running
the workflow.

## Running a workflow

### Prerequisites

To compile the C++ code the requirements for it have to be fulfilled:
1. Cmake (>3.20)
2. PCAP (optional)
3. QRNG library. (optional)
4. Boost.

### Execution

1. To run _a_ workflow first create a directory within which the workflow will be executed.
2. Copy the appropriate snakefile, YAML configuration file and the entire `scripts` directory into the
newly created directory.
3. Rename the copied Snakefile, e.g. `Snakefile_PRNG.smk`, into `Snakefile.smk`.
4. Load the Python virtual environment in which Snakemake has been installed.
5. Execute `snakemake --cores=1 -s Snakefile.smk run/1/single/completed` to run the workflow using a single thread for the first entropy file/seed with
the fundamental data type of the simulation being `float`. To run the workflow for multiple seeds and both precisions execute
`snakemake --cores=1 -s Snakefile.smk run/{1..10}/{single,double}/completed`.

The workflows can be run in parallel by increasing the number of cores. This number dictates how many of the `collectData` rules
will be running in parallel. Note that each assay script will be using the number of threads provided in the configuration file.
It is thus advisable to ensure that the number of cores * number of threads < number of CPU cores in the system.

For the QRNG workflow, whether using `/dev/random` or a dedicated QRNG, it is advisable to use `--cores=1` to prevent contention
for the devices!

## Remarks

In the current state the (p)PRNG and QRNG workflows are essentially duplicates and will have to be modularised and generalised
in the future, to avoid copy-paste errors.

**WARNING**:
The workflow currently copies the _entire repository_ from the base path into the appropriate
subdirectory of 'src'. This is done due to the way CMakeLists.txt is configured to include the compiled
libraries. A run of the current configuration for single & double-precision RNs
will result in the workflow directory growing up to 3GiB in size.
