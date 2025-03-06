/* Monte Carlo Pi Calculation.
 *
 * This file implements an approximation of Pi via
 * Direct Sampling Monte Carlo (DSMC) of the unit disk
 * on the unit square [0,1] x [0,1].
 */

#include <cstdio>
#include <iostream>
#include <iomanip>

#include <chrono>
#include <cstdint>
#include <ios>

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
    // Step 2: Open the target file for writing
    std::ofstream resultFile;
    bool storeMeasurement ( false );
    if (Par.outputFilename != "")
        {
        resultFile.open ( Par.outputFilename, std::ios_base::binary );
        storeMeasurement = true;
        }
    else
        {
        std::cerr << "ERROR! Missing a file name for the output!" << std::endl;
        std::cerr << "Aborting ... " << std::endl;
        return -1;
        }
    // Step 3: Reserve a buffer array to draw the integers into
    // This will be the publicly accessible copy of the QRNGDeviceClasses
    // internal buffer
    QRNGRandomDevice::result_type *outputBuffer = new QRNGRandomDevice::result_type[bufferSize];
    // Step 4: Initialise the QRNG generator
    QRNGRandomDevice generator ( bufferSize );
    // Step 5: Determine how many times the buffer will have to be filled
    // to obtain the desired file size
    std::uint32_t bufferSizeInCharacters = sizeof (QRNGRandomDevice::result_type ) * bufferSize;
    const std::uint32_t charSizeInBytes ( 1 );
    std::uint32_t bufferSizeInBytes = bufferSizeInCharacters * charSizeInBytes;
    std::uint32_t requiredFullBuffers = Par.filesize / bufferSizeInBytes;//Integer division!
    std::uint32_t remainingBytes = Par.filesize - requiredFullBuffers * bufferSizeInBytes;
    
    std::cout << "Using the fundamental datatype: " << fundamentalType.name() << std::endl;
    std::cout << "Size of the buffer type " << sizeof (QRNGRandomDevice::result_type ) << std::endl;
    std::cout << "Requested file size " << Par.filesize << " [bytes] " << std::endl;
    std::cout << "Used buffer size " << bufferSizeInBytes << " [bytes] " << std::endl;
    std::cout << "Required buffer fillings " << requiredFullBuffers << std::endl;
    std::cout << "Byte remnant " << remainingBytes << "[bytes]" << std::endl;
    // Step 6: Fill & store full buffers
    for (std::uint32_t repetition ( 0 );repetition < requiredFullBuffers;++repetition)
        {
        std::cout << "Buffer filling " << repetition + 1 << " of " << requiredFullBuffers << std::endl;
        // fill
        for (uint el ( 0 );el < bufferSize;++el)
            outputBuffer[el] = generator ( );
        //store
        resultFile.write ( reinterpret_cast<char*> ( outputBuffer ), bufferSizeInCharacters );
        //blank the buffer for good measure
        for (uint el ( 0 );el < bufferSize;++el)
            outputBuffer[el] = static_cast<QRNGRandomDevice::result_type> ( 0 );
        }
    // Step 7: Fill and store the remnant, if necessary
    if (remainingBytes != 0)
        {
        uint remainingElements = remainingBytes / sizeof ( QRNGRandomDevice::result_type );
        std::cout << "Remaining elements: " << remainingElements << std::endl;
        for (uint j ( 0 );j < remainingElements;++j)
            outputBuffer[j] = generator ( );
        resultFile.write ( reinterpret_cast<char*> ( outputBuffer ), remainingBytes );
        }
    resultFile.close ( );
    delete[] outputBuffer;

    return 0;
    }
