/* 
 * File:   supportFunctions.h
 * Author: Anton Lebedev
 *
 * Created on 26. April 2023

 */

#pragma once

#include <random>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <typeinfo>

using uint=unsigned int;
#ifdef SINGLE_PRECISION
using dtype=float;
#else
using dtype=double;
#endif

/*! \brief Datatype used to store parameters to be read in.
 * 
 * The structure is used in conjunction with BOOST program_options library
 * to allow the user to vary the parameters dependent on their requirements.
 */
struct params {
    std::string outputFilename;         //!< Name of the file containing the output.
    uint dim;                           //!< Number of dimensions to sample.
    uint numSamplePoints;               //!< Number of sample points of chosen dimension.
    int initSeed;                       //!< Seed for the PRNG, or <0 for time-seeding.
};

/*! \brief Function for parsing the program options
 *
 * \param[in] argc - Number of command-line arguments
 * \param[in] argv - Array containing argument line parameters.
 * \param[in,out] par - Variable which will contain the initialised parameters.
 * 
 * The function takes the command-line arguments via argc, **argv and
 * uses the BOOST program_options library to parse them and store
 * the values into the par variable.
 */
void prog_opt_init(const int argc, const char *const argv[], params& par);

template<typename valueType>
void storeMeasurementsToDisk(std::ofstream& file, const valueType *measurements, const uint Nmeasurements){
    //determine the number of measurements, then store them.
    uint outputDigits=std::numeric_limits<valueType>::digits10;
    //uint Nmeasurements = measurements.size();

    //store number of measurements, mean, and variance
    //file << Nmeasurements << ", " << std::setprecision(outputDigits) << mean << ", "  << std::setprecision(outputDigits) << variance << ", " ;
    //store the separate measurements
    for(uint j(0); j<Nmeasurements; ++j)
        file << std::setprecision(outputDigits) << measurements[j] << ", ";

    file << std::endl;
}
