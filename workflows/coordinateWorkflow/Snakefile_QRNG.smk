configfile: "coord_QRNG_snakemake.yaml"

precisions=['single', 'double']
entropyFiles=[idx for idx in range(1,10+1)]

def getEntropyFiles(wildcards):
    return config["entropy"][wildcards.entropyFileIndex]

def getRepoPath(wildcards):
    print(config['repoPath'])
    return config['repoPath']

def getCompileFlags(wildcards):
    if wildcards.precision=='single':
        return "-DSINGLE_PRECISION"
    else:
        return " "

# compilation
rule copyRepo:
    input:
        srcDir=getRepoPath
    message: "Copying source from {input.srcDir}"
    output:
        directory("src/{entropyFileIndex}")
    shell:
        "cp -R {input.srcDir} {output}"

rule runCompilation:
    input:
        rules.copyRepo.output,
    message: "{wildcards.precision} {wildcards.entropyFileIndex}"
    params:
        precisionOption=getCompileFlags,
        entropyFile=getEntropyFiles
    output:
        directory("build/{entropyFileIndex}/{precision}")
    shell:
        "sed -i 's|/dev/random|{params.entropyFile}|g' {input}/DSMCPI/include/QRNGRandomDevice.h && "
        "mkdir {output} && "
        "cd {output} && "
        "cmake -DCMAKE_CXX_FLAGS='{params.precisionOption} -DQRNG -DQRNG_USE_DEVICE=0' ../../../{input}/nDCoordinates &&"
        "make"

rule collectBinaries:
    input:
        "build/{entropyFileIndex}/{precision}"
    message: "{wildcards.precision} {wildcards.entropyFileIndex}"
    output:
        "bin/{entropyFileIndex}/{precision}/coordinates"
    shell:
        "cp {input}/coordinates {output}"

# Execution

rule collectDataPoints:
    input:
        "bin/{entropyFileIndex}/{precision}/coordinates"
    output:
        "run/{entropyFileIndex}/{precision}/coordinates.lst"
    params:
        numSamples=config['samples'],
        dim=config['dim']
    threads: 1
    shell:
        "./{input} --numSamplePoints {params.numSamples} --dimension {params.dim} --outfile {output}"

# visualisation & batch processing

rule visualiseNNDistributions:
    input:
        "run/{entropyFileIndex}/{precision}/coordinates.lst"
    params:
        groups=config['batches'],
        minY=config['minY'],
        maxY=config['maxY']
    output:
        "vis/nndist_{precision}_{entropyFileIndex}.png",
        "vis/nndata_{precision}_{entropyFileIndex}.dat",
        "vis/nnStatistics_{precision}_{entropyFileIndex}.dat",
        "vis/nndata_PE_Q_{precision}_{entropyFileIndex}.dat",
        "vis/nndata_Rm_{precision}_{entropyFileIndex}.dat",
        "vis/nndata_Avg_{precision}_{entropyFileIndex}.dat",
        "vis/nndata_CoV_{precision}_{entropyFileIndex}.dat",
    script:
        "scripts/nearestNeighbourPlot.py"

rule visualiseLESDistributions:
    input:
        "run/{entropyFileIndex}/{precision}/coordinates.lst"
    params:
        groups=config['batches'],
        minY=config['minY'],
        maxY=config['maxY']
    output:
        "vis/lesdist_{precision}_{entropyFileIndex}.png",
        "vis/lesdata_{precision}_{entropyFileIndex}.dat",
        "vis/lesStatistics_{precision}_{entropyFileIndex}.dat",
    script:
        "scripts/LESplot.py"

# plot statistics of the collected means and variances
# aggregate all datasets into one plot.

rule plotMeanLESStatistics:
    input:
        expand("vis/lesStatistics_{precision}_{entropyFileIndex}.dat", precision=precisions, entropyFileIndex=entropyFiles )
    output:
        "vis/lesMeanValueStatistics.png"
    script:
        "scripts/meanDistribution.py"

rule plotMeanNNStatistics:
    input:
        expand("vis/nnStatistics_{precision}_{entropyFileIndex}.dat", precision=precisions, entropyFileIndex=entropyFiles )
    output:
        "vis/nnMeanValueStatistics.png"
    script:
        "scripts/meanDistribution.py"
