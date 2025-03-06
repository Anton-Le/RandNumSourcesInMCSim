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
    const uint dim ( Par.dim ), nSamples ( Par.numSamplePoints );

    //Step 4: Reserve a plain array to store coordinates to
    dtype **coordinates = new dtype*[nSamples];
    for(uint j(0); j<nSamples; ++j){
        coordinates[j] = new dtype[dim];
        for(uint i(0); i<dim; ++i)
            coordinates[j][i] = static_cast<dtype>(0);
    }

    std::cout << "Running coordinate collection for " << nSamples << " vectors" << std::endl;
    std::cout << "of dimension " << dim << std::endl;
    std::cout << "Each coordinate is sampled from [0,1]" << std::endl;
    std::cout << "For a total of " << nSamples * dim << " random numbers." << std::endl;
    std::cout << "Using the fundamental datatype: " << fundamentalType.name() << std::endl;
    //Step 5: Open a file for data storage, if requested
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
    //Step 6: Loop over all samples in a parallelised fashion
#pragma omp parallel shared(coordinates, rngSeed)
    {
        rngInitStart = std::chrono::system_clock::now ( );
        int threadId = omp_get_thread_num();
#ifndef QRNG
        std::cout << "Using PRNG" << std::endl;
        std::mt19937 generator;
        //std::ranlux48 generator;
        generator.seed(rngSeed + threadId);
#else
        //-------------- QRNG - INIT
        std::cout << "Using QRNG" << std::endl;
        QRNGRandomDevice generator(1024*1024);
#endif
        std::uniform_real_distribution<dtype> prndist ( 0.0, 1.0 );
        auto randomNumberDraw = std::bind ( prndist, std::ref ( generator ) );
        rngInitStop = std::chrono::system_clock::now ( );


        //Step 6.1: Loop over all coordinates 
#pragma omp for collapse(2)
        for( uint sampleId = 0; sampleId < nSamples; ++sampleId){
            for( uint d = 0; d < dim; ++d){
                coordinates[sampleId][d] = randomNumberDraw();
                }
            }
#pragma omp barrier 
    }
    rngInitDuration = rngInitStop - rngInitStart;
    batchLoopStop = std::chrono::system_clock::now ( );
    auto batchLoopDuration = batchLoopStop - batchLoopStart;
    //Step78: Store data to disk
    storageStart = std::chrono::system_clock::now ( );
    if (storeMeasurement){
        for (uint sampleId ( 0 );sampleId < nSamples;++sampleId)
            {
            //store the estimates
            storeMeasurementsToDisk ( resultFile, coordinates[sampleId], dim );
            }
        }
    //Step 8: Close file if opened
    if (storeMeasurement)
        resultFile.close ( );
    storageStop = std::chrono::system_clock::now ( );
    storageDuration = storageStop - storageStart;
    totalStorageDuration = std::chrono::duration_cast<std::chrono::microseconds>( storageDuration ).count ( );
    
    //Step 9: clean-up of the measurement data
    for (uint sampleId ( 0 );sampleId < nSamples;++sampleId)
        {
        delete[] coordinates[sampleId];
        }
    delete[] coordinates;

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
