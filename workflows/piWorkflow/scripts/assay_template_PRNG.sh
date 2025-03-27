#!/bin/bash

#
# Run a tensor-product parameter space check
# Assuming the revision of the code which runs batches of simulations
# Using: Nsamples = [1000, 10000, 100000, 1000000] samples per Pi/4 calculation
# Running: Nrep = [2,8,32,64,256,1024] repetitions of each Pi/4 calculation to obtain the mean
# Repeating the procedure: Nbatches = [2, 8, 16, 64, 128] timesto obtain the mean of means
#
ORIGINDIR=$(pwd)
PROGNAME=EXECNAME
RESULT_FILE_PREFIX="results_assay_ASSAYID_"


#Assay parameters
Nsamples=(SAMPLES_PER_REP)
Nrep=(REP_PER_BATCH)
Nbatches=(BATCHES)
Method=METHOD
Seed=SEED

export OMP_NUM_THREADS=threads

echo "Source directory: "$ORIGINDIR
echo "Program directory: "$PROGDIR

# Print assay parameters
echo "Running assay with a fixed number of points per repetition"
echo "Total number of samples to approximate pi: "$NptsTotal
echo "Repetitions per batch: "${Nrep[@]}
echo "Points per repetition: "$NptsPerRepetition

#-----------------------------------------
cp $PROGDIR/$PROGNAME $ORIGINDIR/$PROGNAME
echo "running pi approximation"
for batchSize in "${Nbatches[@]}"; do
    echo "Running with "$batchSize" batch size";
        for rep in "${Nrep[@]}"; do
            echo "Running batch with "$rep" number of repetitions"
            for samples in "${Nsamples[@]}"; do
                RESULT_FILE_NAME=$RESULT_FILE_PREFIX"_batchSize_"$batchSize"_Nrep_"$rep"_samples_"$samples".txt"
                $PROGNAME --initSeed $Seed --method $Method --nBatches $batchSize --nRep $rep --samplesPerRepetition $samples --outfile $RESULT_FILE_NAME
	            wait
            done
            wait
        done
	    wait
done
