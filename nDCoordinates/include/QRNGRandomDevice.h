/*
 * The MIT License
 *
 * Copyright 2023 Anton Lebedev.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* 
 * File:   QRNGRandomDevice.h
 * Author: Anton Lebedev
 *
 * Created on 18. Mai 2023, 12:03
 */



#ifndef QRNGRANDOMDEVICE_H
#define QRNGRANDOMDEVICE_H

#include <cstdint>
#include <fstream>
#include <byteswap.h>

#if QRNG_USE_DEVICE == 1
#include "qrng_api.h"
#else
/**
 * Proxy status code enum to account for the absence of a QRNG device.
  */
  
typedef enum {
	QRNG_ERROR_NET_RETRIES_EXCEEDED = -14,
	QRNG_ERROR_INSUFFICIENT_ENTHROPY = -10,
	QRNG_ERROR_READING_DEVICE = -7,
	QRNG_ERROR_INCOMPLETE_DATA = -6,
	QRNG_ERROR_OPENING_DEVICE = -5,
	QRNG_NOT_INITIALIZED = -4,
	QRNG_ALREADY_INITIALIZED = -3,
	QRNG_NO_DEVICE_FOUND = -1,
	QRNG_SUCCESS = 0,
}Qrng_status;
// Typedefs adapted to conform with the QRNG api
typedef signed char s8;
typedef unsigned char u8;
typedef signed short s16;
typedef unsigned short u16;
typedef signed long s32;
typedef unsigned long u32;
typedef signed long long s64;
typedef unsigned long long u64;
#endif


namespace qrng_random{
        std::ifstream datastream;
}


/*! \brief A proxy function for qrng_init.
 * 
 * This function opens a file stream from
 * `/dev/urandom` to be used in subsequent functions.
 *
 * \return Error code
 */
int qrng_init_proxy(){
    if( qrng_random::datastream.is_open() == true)
        return -3;

    // Warn the user that proxy functions are being used
    std::cerr << "WARNING!!! Using file-based proxy for the entropy source!" << std::endl;

    qrng_random::datastream.open("/dev/random", std::fstream::in | std::fstream::binary);
    if( qrng_random::datastream.is_open() != true)
        return -5;

    return 0;
}

/*! \brief A proxy function for qrng_get_status.
 * 
 * This function checks whether any error flags were
 * set for the data stream
 *
 * \return Error code
*/
int qrng_get_status_proxy(){
    if( qrng_random::datastream.good() == false )
        return -7;
    
    //good behaviour:
    return 0;
}

/*! \brief A proxy function for `qrng_deinit`.
 *
 * This function closes the datastream if it's open,
 * otherwise it returns an error code similar to `qrng_deinit`.
 *
 * \return Error code of the operation (0 on success).
 */
int qrng_deinit_proxy(){
    if( qrng_random::datastream.is_open() != true)
        return -4;

    qrng_random::datastream.close();

    return 0;
}

/*! \brief A proxy function for qrng_get.
 *
 * This function posesses the same signature as `qrng_get`
 * but uses a stream from `/dev/urandom` as entropy source.
 * \param[in,out] data - array to store random bits to.
 * \param[in] requestedCharacters - number of bytes requested.
 * \param[out] obtainedCharacters - number of bytes actually obtained.
 * \return Default value of 0.
 */
int qrng_get_proxy(u8* data, const u32 requestedCharacters, u32 *providedCharacters){
    //allocate an internal buffer
    char* rawData = new char[requestedCharacters];
    qrng_random::datastream.read(rawData, requestedCharacters);
    int obtainedCharacters = qrng_random::datastream.gcount();
    for(uint j(0); j< obtainedCharacters; ++j)
        data[j] = rawData[j];
    qrng_random::datastream.clear();
    delete[] rawData;
    *providedCharacters = obtainedCharacters;
    return 0;
}

/*! \brief A proxy function for entropy counting.
 * 
 * This function returns the entropy of `/dev/random`.
 * It is intended to be used with the default proxy RN
 * source on Linux.
 * **CAVEAT**: Usage of a different binary file as an
 * entropy source or use of a different OS will result
 * in indeterminate behaviour.
 */
int qrng_get_entropy(){
    std::ifstream entropyCountFile;
    entropyCountFile.open("/proc/sys/kernel/random/entropy_avail", std::fstream::in);
    if( entropyCountFile.is_open() != true)
        return 0;
    else{
        int entropyAvail;
        entropyCountFile >> entropyAvail;
        return entropyAvail;
    }
}


/**
 * A standard interface to a buffer-type entropy source.
 * 
 * The class contains an array of unsigend integers of default size 10 MiB
 * as an entropy buffer and a private `fill` function which
 * calls `qrng_get` within an OpenMP `critical` region to ensure
 * that only one OpenMP thread at a time can fetch numbers from the device.
 * 
 */
class QRNGRandomDevice {
public:
    /** The type of the generated random value. */
    using result_type = unsigned int;
    using counterType = unsigned int;

    // constructors, destructors and member functions
    QRNGRandomDevice(const counterType bufferSize = 10485760 ) : \
        bufferSizeInNumbers(bufferSize), availableBits(0), \
        qrngState(Qrng_status::QRNG_SUCCESS), rnCounter(0) {
        this->requiredBytes = sizeof (result_type);
        this->entropyBuffer = new result_type[this->bufferSizeInNumbers];

#if QRNG_USE_DEVICE == 1
        // temporary hard-coded init parameters
        Qrng_init_param deviceInitParameters = { Qrng_board_type::QRNG_VERTEX_A1,\
                Qrng_interface_type::QRNG_10G_NET, static_cast<char*>("enp3s0") };
        //old device
        //Qrng_init_param deviceInitParameters = { Qrng_board_type::QRNG_APEX_STREAM,\
	    //   	Qrng_interface_type::QRNG_PCIE };
#endif
        //initialise the QRNG if it hasn't already happened.
        //There's no harm in calling `qrng_init` multiple times in series.
#pragma omp critical
        {

#if QRNG_USE_DEVICE == 1
            this->qrngState = qrng_init_param(deviceInitParameters) ;
#else
            this->qrngState = qrng_init_proxy();
#endif

        if( (this->qrngState != Qrng_status::QRNG_SUCCESS) \
            && (this->qrngState != Qrng_status::QRNG_ALREADY_INITIALIZED) \
	    && ( this->qrngState != Qrng_status::QRNG_ERROR_NET_RETRIES_EXCEEDED ) ){
            std::cerr << "Error initialising the device" << std::endl;
            std::cerr << "Error code: " << this->qrngState << std::endl;
            exit(1);
        }
        }
    };
    /*! \brief Destructor.
     * 
     * The destructor of the class frees the entropy buffer for each object.
     * Only if the object is part of the `master` OpenMP thread will `qrng_deinit`
     * be called, to ensure it's called exactly once.
     * 
     * This implementation **assumes** that the threads are synchronized before
     * leaving the parallel region, else the QRNG state will be indeterminate!
     */
    ~QRNGRandomDevice() {
        delete[] this->entropyBuffer;
#pragma omp master
        {
#if QRNG_USE_DEVICE == 1
            int qrngState = qrng_get_status();
#else
            int qrngState = qrng_get_status_proxy();
#endif
            if(qrngState != Qrng_status::QRNG_NOT_INITIALIZED)
#if QRNG_USE_DEVICE == 1
                this->qrngState = qrng_deinit();
#else
                this->qrngState = qrng_deinit_proxy();
#endif
        }
    }

    static constexpr result_type
    min() {
        return std::numeric_limits<result_type>::min();
    }

    static constexpr result_type
    max() {
        return std::numeric_limits<result_type>::max();
    }

    double
    entropy() const noexcept {
        //Query the system
        
#if QRNG_USE_DEVICE == 1
        return 0.0;
#else
        return qrng_get_entropy();
#endif
    }

    result_type
    operator()() {
        while (this->availableBits < 1) {
            //refill the buffer when it's been emptied
            this->availableBits = this->fillBuffer(this->bufferSizeInNumbers);
            this->rnCounter = 0;
        }
        //decrement the number of available RNs
        this->availableBits--;
        return entropyBuffer[this->rnCounter++];

    }

    // No copy functions.
    QRNGRandomDevice(const QRNGRandomDevice&) = delete;
    void operator=(const QRNGRandomDevice&) = delete;

private:
    /*! \brief Function to fill the entropy buffer from the device.
     * 
     * This function determines how many number of the type `qrng_get`
     * provides correspond to the `requestedNumbers` amount of
     * the entropy buffer type and will request as many.
     * Upon completion the function returns the number of
     * RNs of `result_type` actually obtained.
     */
    counterType fillBuffer(const counterType requestedNumbers) {
        //std::cout << " Filling the buffer with RNs... Requested: " << requestedNumbers << std::endl;
        std::uint8_t QRNGdataTypeSize = sizeof(u8);
        std::uint8_t returnDataTypeSize = sizeof(result_type);
        std::uint32_t requestedRNsOfQRNGType \
            = requestedNumbers * returnDataTypeSize / QRNGdataTypeSize;
        std::uint32_t obtainedRNsOfDataType(0);
        u32 providedRNsOfQRNGType(0);

#if QRNG_USE_DEVICE == 1
#pragma omp critical
        this->qrngState = qrng_get(reinterpret_cast<u8*>(entropyBuffer), \
                        static_cast<u32>(requestedRNsOfQRNGType),
                        &providedRNsOfQRNGType);

#else
#pragma omp critical
        this->qrngState = qrng_get_proxy(reinterpret_cast<u8*>(entropyBuffer), \
                        static_cast<u32>(requestedRNsOfQRNGType),
                        &providedRNsOfQRNGType);
#endif
            
        if (this->qrngState != Qrng_status::QRNG_SUCCESS) {
                std::cerr << "[QRNG] Error drawing numbers!" << std::endl;
                std::cerr << "[QRNG] Error code: " << this->qrngState << std::endl;
                //it's probably better to throw an exception here
                //exit(1);
        }
        obtainedRNsOfDataType = static_cast<counterType>(providedRNsOfQRNGType * QRNGdataTypeSize / returnDataTypeSize);
        return obtainedRNsOfDataType;
    }

    counterType bufferSizeInNumbers;
    counterType availableBits;
    counterType rnCounter;
    counterType requiredBytes;
    int qrngState;

    result_type *entropyBuffer;

};

#endif /* QRNGRANDOMDEVICE_H */

