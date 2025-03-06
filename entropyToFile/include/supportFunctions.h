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

//!< Buffer size in integers
constexpr uint bufferSize(8000000);

/*! \brief Datatype used to store parameters to be read in.
 * 
 * The structure is used in conjunction with BOOST program_options library
 * to allow the user to vary the parameters dependent on their requirements.
 */
struct params {
    std::string outputFilename;         //!< Name of the file containing the output.
    uint filesize;                      //!< File size in Bytes
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
