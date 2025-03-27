/* Monte Carlo Pi Calculation.
 *
 * This file implements an approximation of Pi via
 * Direct Sampling Monte Carlo (DSMC) of the unit disk
 * on the unit square [0,1] x [0,1].
 */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <functional>
#include <cstdint>
#include <omp.h>
#include <trng/mrg5s.hpp>

#include "supportFunctions.h"
#include "QRNGRandomDevice.h"

int
main ( int argc, char** argv )
    {
    //Step 0: Store the type information on the base datatype
    const std::type_info& fundamentalType = typeid(dtype);
    //Step 1: Load program parameters
    params Par;
    prog_opt_init ( argc, argv, Par );

    //Step 2: Initialise the RNG in accordance with prescription
    uint rngSeed ( 0 );
    if (Par.initSeed < 0.0)
        {
        std::chrono::time_point<std::chrono::system_clock> initTime ( std::chrono::system_clock::now ( ) );
        rngSeed = std::chrono::duration_cast<std::chrono::microseconds>( initTime.time_since_epoch ( ) ).count ( );
        }
    else
        {
        rngSeed = Par.initSeed;
        }

    std::chrono::time_point< std::chrono::system_clock > rngInitStart, rngInitStop;

    auto rngInitDuration = rngInitStop - rngInitStart;

    //Step 3: Retrieve frequently used parameters for local use.
    const uint nBatches ( Par.nBatches ), nRep ( Par.nRep );
    const uint nSamplesPerRepetition ( Par.samplesPerRepetition );



    //Step 5.1: Reserve an array for all pi samples
    dtype **batchSampleMeans = new dtype*[nBatches];
    dtype **batchSampleAbsError = new dtype*[nBatches];
    std::vector<dtype> batchMean ( nBatches, 0.0 );
    std::vector<dtype> batchVariance ( nBatches, 0.0 );
    std::vector<dtype> batchMeanAbsErr ( nBatches, 0.0 );
    std::vector<dtype> batchAbsErrVariance ( nBatches, 0.0 );

    for (uint j ( 0 );j < nBatches;++j)
        {
        batchSampleMeans[j] = new dtype[nRep];
        batchSampleAbsError[j] = new dtype[nRep];

        for (uint i ( 0 );i < nRep;++i)
            batchSampleMeans[j][i] = batchSampleAbsError[j][i] = 0.0;
        }


    std::cout << "Running data collection for " << nBatches << " batches" << std::endl;
    std::cout << "of " << nRep << " pi-estimations per batch" << std::endl;
    std::cout << "each using " << nSamplesPerRepetition << " sample points on [0,1] x [0,1]" << std::endl;
    std::cout << "For a total of " << nSamplesPerRepetition * nRep * nBatches << " sample points in [0,1] x [0,1]" << std::endl;
    std::cout << "Using the fundamental datatype " << fundamentalType.name() << std::endl;
    //Step 6: Open a file for data storage, if requested
    std::ofstream resultFile;
    bool storeMeasurement ( false );
    if (Par.outputFilename != "")
        {
        resultFile.open ( Par.outputFilename );
        storeMeasurement = true;
        }

    std::chrono::time_point< std::chrono::system_clock > batchLoopStart, batchLoopStop;
    std::chrono::time_point< std::chrono::system_clock > storageStart, storageStop;
    auto storageDuration = storageStop - storageStart;
    std::uint32_t totalStorageDuration ( 0 );
    batchLoopStart = std::chrono::system_clock::now ( );
    //Step 7: Loop over all batches
#pragma omp parallel shared(batchSampleMeans, batchSampleAbsError, rngSeed, batchMean, batchVariance, batchMeanAbsErr, batchAbsErrVariance)
    {
        rngInitStart = std::chrono::system_clock::now ( );
        const uint numThreads(omp_get_num_threads());
        uint threadId( omp_get_thread_num() );
#ifndef QRNG
        std::cout << "Using PRNG" << std::endl;
        trng::mrg5s generator;
        generator.seed(static_cast<unsigned long>(rngSeed) );
        // sequence splitting
        // generator.split(numThreads, threadId);
        // block splitting
        generator.jump( 2 * (threadId * nSamplesPerRepetition * nRep * nBatches / numThreads) );
#else
        //-------------- QRNG - INIT
        std::cout << "Using QRNG" << std::endl;
        QRNGRandomDevice generator;
#endif
        
        std::uniform_real_distribution<dtype> prndist ( 0.0, 1.0 );
        if( Par.method == 2)
            prndist.param( std::uniform_real_distribution<dtype>::param_type(0.0, 0.5) );
        
        auto randomNumberDraw = std::bind ( prndist, std::ref ( generator ) );
        rngInitStop = std::chrono::system_clock::now ( );
        //Step 4: reserve space for the markers
        std::vector<bool> samples ( nSamplesPerRepetition );
        std::fill ( samples.begin ( ), samples.end ( ), false );

        //Step 5: Reserve an array for sample means and variances
        std::vector<dtype> sampleMeans ( nRep ), absError ( nRep );
        std::fill ( sampleMeans.begin ( ), sampleMeans.end ( ), 0.0 );
        std::fill ( absError.begin ( ), absError.end ( ), 0.0 );
        dtype exptMean ( 0 ), exptVariance ( 0 );
        uint withinUnitCircle ( 0 );
#pragma omp for 
        for (uint batchId = 0;batchId < nBatches;++batchId)
            {
            //Step 7.1: Perform nRep repetitions of the Pi approximation
            for (uint repetition = 0;repetition < nRep;++repetition)
                {
                //Step 7.1.1: Draw the samples and check condition
                if(Par.method == 2){
                    sampleAndCheckPiCondition ( samples, nSamplesPerRepetition, randomNumberDraw );
                }
                else{
                    sampleUnitSquareAndCheckPiCondition ( samples, nSamplesPerRepetition, randomNumberDraw );
                }
                //Step 7.1.2: Compute the approximation via counter
                withinUnitCircle = std::count ( samples.begin ( ), samples.end ( ), true );
                if( Par.method == 2 ){
                    sampleMeans[repetition] = 2.0 / (static_cast<dtype> ( withinUnitCircle ) / static_cast<dtype> ( nSamplesPerRepetition ) );
                }
                else{
                    sampleMeans[repetition] = 4.0 * static_cast<dtype> ( withinUnitCircle ) / static_cast<dtype> ( nSamplesPerRepetition );
                }
                absError[repetition] = std::fabs ( sampleMeans[repetition] - pi );
                //store to global array to be stored to file
                batchSampleMeans[batchId][repetition] = sampleMeans[repetition];
                batchSampleAbsError[batchId][repetition] = absError[repetition];
                }
            //Step 7.2: Determine the mean value of the Pi estimates and their abs. errors
            exptMean = sampleMean ( sampleMeans );
            dtype meanAbsErr = sampleMean ( absError );
            //Step 7.3: Determine the empirical variance
            exptVariance = sampleVariance ( sampleMeans, exptMean );
            dtype absErrVariance = sampleVariance ( absError, meanAbsErr );
            //Store data for en-block storage
            batchMean[batchId] = exptMean;
            batchVariance[batchId] = exptVariance;
            batchMeanAbsErr[batchId] = meanAbsErr;
            batchAbsErrVariance[batchId] = absErrVariance;
            //Step 7.4: Print data for the batch
            std::cout << "RESULTS of batch " << batchId << std::endl;
            std::cout << "===============================" << std::endl;
            std::cout << "N = " << nBatches << std::endl;
            std::cout << "<pi> = " << std::setprecision ( 16 ) << exptMean << std::endl;
            std::cout << "Var[<pi>] = " << std::setprecision ( 16 ) << exptVariance << std::endl;
            std::cout << "< |pi_j - pi| > = " << std::setprecision ( 16 ) << meanAbsErr << std::endl;
            std::cout << "Var[ |pi_j - pi| ] = " << std::setprecision ( 16 ) << absErrVariance << std::endl;
            }
        #pragma omp barrier
        rngInitDuration = rngInitStop - rngInitStart;
    }
    batchLoopStop = std::chrono::system_clock::now ( );
    auto batchLoopDuration = batchLoopStop - batchLoopStart;
    //Step 8: Store data to disk
    storageStart = std::chrono::system_clock::now ( );
    if (storeMeasurement)
        for (uint batchId ( 0 );batchId < nBatches;++batchId)
            {
            //store the estimates
            storeMeasurementsToDisk ( resultFile, batchMean[batchId], batchVariance[batchId], batchSampleMeans[batchId], nRep );
            //store the absolut errors
            storeMeasurementsToDisk ( resultFile, batchMeanAbsErr[batchId], batchAbsErrVariance[batchId], batchSampleAbsError[batchId], nRep );
            }
    storageStop = std::chrono::system_clock::now ( );
    storageDuration = storageStop - storageStart;
    totalStorageDuration = std::chrono::duration_cast<std::chrono::microseconds>( storageDuration ).count ( );
    //clean-up of the measurement data
    for (uint j ( 0 );j < nBatches;++j)
        {
        delete[] batchSampleMeans[j];
        delete[] batchSampleAbsError[j];
        }
    delete[] batchSampleAbsError;
    delete[] batchSampleMeans;
    //Step 8: Close file if opened
    if (storeMeasurement)
        resultFile.close ( );
    //Step 9: Print the timings
    std::cout << "Time to initialise the RNG: " \
            << std::chrono::duration_cast<std::chrono::microseconds>( rngInitDuration ).count ( ) \
            << " [µs]" << std::endl;
    std::cout << "Total simulation time: " \
            << std::chrono::duration_cast<std::chrono::microseconds>( batchLoopDuration ).count ( ) \
            << " [µs]" << std::endl;
    std::cout << "Of which " << totalStorageDuration << " [µs] were spent storing to file." << std::endl;
    return 0;
    }
