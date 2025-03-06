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

constexpr dtype pi(3.1415926535897932384626433);
/*! \brief Datatype used to store parameters to be read in.
 * 
 * The structure is used in conjunction with BOOST program_options library
 * to allow the user to vary the parameters dependent on their requirements.
 */
struct params {
    std::string outputFilename;         //!< Name of the file containing the output.
    uint nBatches;                      //!< Number of batches to simulate.
    uint nRep;                          //!< Number of repetitions per batch.
    uint samplesPerRepetition;          //!< Number of samples to draw per repetition.
    int initSeed;                       //!< Seed for the PRNG, or <0 for time-seeding.
    uint method;                        //!< Selection of the method (standard or Buffon)
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

/*! \brief Function to sample the unit square and check the Pi-condition.
 *
 * \param[in,out] markers - vector of booleans to store the results of the checks in.
 * \param[in] nSamples - number of samples to be distributed on [0,1]^2.
 * \param[in,out] rng - random-number generating function.
 * 
 * It is assumed that the rng argument carries an internal state and is hence
 * formally an in/out parameter. The function utilises it only as an incoming parameter.
 */
template <typename Func>
void sampleUnitSquareAndCheckPiCondition(std::vector<bool>& markers, const uint nSamples, Func& rng ){
	std::vector<dtype> x(nSamples, static_cast<dtype>(0)), y(nSamples, static_cast<dtype>(0));
	// fill the coordinates
	for(uint sampleId(0); sampleId < nSamples; ++sampleId){
	   x[sampleId] = rng();
	//} //if these two lines are commented out we perform XYXY sampling, otherwise XXYY sampling
	//for(uint sampleId(0); sampleId < nSamples; ++sampleId){
	   y[sampleId] = rng();
	}
	// clear markers
        markers.assign( nSamples, false );
	for(uint sampleId(0); sampleId < nSamples; ++sampleId){
	   if( x[sampleId]*x[sampleId]+y[sampleId]*y[sampleId] <= static_cast<dtype>(1.0) )
	       markers[sampleId] = true;
	}
}

/*! \brief Function to sample [0,1/2] x [0, Pi/2] and check for intersections.
 *
 * \param[in,out] markers - vector of booleans to store the results of the checks in.
 * \param[in] nSamples - number of samples to be distributed on [0,1/2] x [0, pi/2] .
 * \param[in,out] rng - random-number generating function, generating RNs uniformly on [0,1/2].
 * 
 * It is assumed that the rng argument carries an internal state and is hence
 * formally an in/out parameter. The function utilises it only as an incoming parameter.
 * The function simulates Buffon's problem to determine a stochastic approximation
 * to 2/Pi.
 */
template <typename Func>
void sampleAndCheckPiCondition(std::vector<bool>& markers, const uint nSamples, Func& rng ){
	std::vector<dtype> x(nSamples, static_cast<dtype>(0)), alpha(nSamples, static_cast<dtype>(0));
	// fill the coordinates
	for(uint sampleId(0); sampleId < nSamples; ++sampleId){
	   x[sampleId] = rng();
	   alpha[sampleId] = pi * rng();
	}
	// clear markers
    markers.assign( nSamples, false );
	for(uint sampleId(0); sampleId < nSamples; ++sampleId){
	   if( dtype(0) <= x[sampleId] && x[sampleId] <= dtype(0.5) * sin(alpha[sampleId]) ) 
	       markers[sampleId] = true;
	}
}

template<typename valueType>
valueType sampleMean(const std::vector<valueType>& samples){
    //Kahan-Neumaier summation
    valueType sum(0), compensation(0), tmp(0);
    for(uint id(0); id < samples.size(); ++id){
        tmp = sum + samples[id];
        if( std::abs( sum ) >= std::abs(samples[id]) ){
            compensation += (sum - tmp) + samples[id];
        }else{
            compensation += (samples[id] - tmp) + sum;
        }
        sum = tmp;
    }
    return (sum+compensation)/static_cast<valueType>(samples.size());
}

template<typename valueType>
valueType sampleVariance(const std::vector<valueType>& samples, const valueType mean){
    std::vector<valueType> deviations(samples.size());
    valueType tmp(0);
    for(uint i(0); i<samples.size(); ++i)
        deviations[i] = (samples[i] - mean) * (samples[i] - mean);
    
    valueType variance = sampleMean(deviations)*static_cast<valueType>(samples.size())/static_cast<valueType>(samples.size() - 1);
    return variance;
}

template<typename valueType>
void storeMeasurementsToDisk(std::ofstream& file, const valueType mean, const valueType variance, const valueType *measurements, const uint Nmeasurements){
    //determine the number of measurements, then store them.
    uint outputDigits=std::numeric_limits<valueType>::digits10;
    //uint Nmeasurements = measurements.size();

    //store number of measurements, mean, and variance
    file << Nmeasurements << ", " << std::setprecision(outputDigits) << mean << ", "  << std::setprecision(outputDigits) << variance << ", " ;
    //store the separate measurements
    for(uint j(0); j<Nmeasurements; ++j)
        file << std::setprecision(outputDigits) << measurements[j] << ", ";

    file << std::endl;
}
