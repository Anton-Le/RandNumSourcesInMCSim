configfile: "configuration.yaml"

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
        "cp -R {input.srcDir}/ToyExamples {output}"


#ruleorder: substituteEntropyFile > copyRepo

rule runCompilation:
    input:
        rules.copyRepo.output,
        #"tmp/{entropyFileIndex}.complete"
    message: "{wildcards.precision} {wildcards.entropyFileIndex}"
    params:
        precisionOption=getCompileFlags,
        entropyFile=getEntropyFiles
    output:
        directory("build/{entropyFileIndex}/{precision}")
    shell:
        "mkdir {output} && "
        "cd {output} && "
        "cmake -DCMAKE_CXX_FLAGS='{params.precisionOption} -march=native' ../../../{input} &&"
        "make"

rule collectBinaries:
    input:
        "build/{entropyFileIndex}/{precision}"
    message: "{wildcards.precision} {wildcards.entropyFileIndex}"
    output:
        "bin/{entropyFileIndex}/{precision}/DSMCPI"
    shell:
        "cp {input}/DSMCPI {output}"
        
######################################
# Execution
######################################

rule copyAssayFiles:
    input:
        "bin/{entropyFileIndex}/{precision}/DSMCPI"
    output:
        "run/{entropyFileIndex}/{precision}/assay.sh"
    shell:
        "cp ./scripts/assay_template_PRNG.sh {output}"

rule instantiateAssayScript:
    input:
        "run/{entropyFileIndex}/{precision}/assay.sh"
    output:
        touch("run/{entropyFileIndex}/{precision}/instantiated")
    params:
        batches=config['batches'],
        samples=config['samples'],
        repetitions=config['repetitions'],
        threads=config['threads'],
        method=config['method'],
        rngSeed=getEntropyFiles
    message: "Running {params.batches} batches with {params.repetitions} per batch and {params.samples} per rep. Using {params.threads} threads"
    shell:
        "sed -i 's|BATCHES|{params.batches}|g' {input} && "
        "sed -i 's|REP_PER_BATCH|{params.repetitions}|g' {input} && "
        "sed -i 's|SAMPLES_PER_REP|{params.samples}|g' {input} && "
        "sed -i 's|threads|{params.threads}|g' {input} &&"
        "sed -i 's|ASSAYID|{wildcards.entropyFileIndex}|g' {input} &&"
        "sed -i 's|EXECNAME|../../../bin/{wildcards.entropyFileIndex}/{wildcards.precision}/DSMCPI|g' {input} && "
        "sed -i 's|METHOD|{params.method}|g' {input} && "
        "sed -i 's|SEED|{params.rngSeed}|g' {input}"

rule collectData:
    input:
        "run/{entropyFileIndex}/{precision}/assay.sh",
        "run/{entropyFileIndex}/{precision}/instantiated"
    output:
        "run/{entropyFileIndex}/{precision}/results_assay_{entropyFileIndex}_*.txt", 
        touch("run/{entropyFileIndex}/{precision}/completed")
    threads: 1
    shell:
        "cd run/{wildcards.entropyFileIndex}/{wildcards.precision} && "
        "./assay.sh"

# visualisation & batch processing