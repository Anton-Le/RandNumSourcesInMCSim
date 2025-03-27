
Exploration of effects of entropy (random numbers) sources on Monte Carlo simulations.

# Direct Sampling Monte Carlo estimation of Pi (DSMCPI)

The `pi-estimation` branch of this repository contains the source
code, [Snakemake](https://snakemake.readthedocs.io/en/stable) workflows and Jupyter notebooks
used to test the effects of different random number generators (RNGs) on simple
Monte Carlo applications.

This branch also contains (c.f. below)
- a tool to dump random bits from the QRNG device (or `/dev/random`) to a binary file and
- a tool to convert those binary files into lists of n-dimensional vectors in $[0,1]^d$.

The latter two (`entropyToFile` and `nDCoordinates`) were written with the intent of 
explaining observations of the DSMCPI estimation 
by allowing an analysis of the actual entropy as well as of the results obtained with it.

### Compilation
In the root directory of the repository issue the following
`mkdir builddir && cd builddir`
`cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-DQRNG -DQRNG_USE_DEVICE=1 -DSINGLE_PRECISION -march=native" ../DSMCPI`
`make`

The flags in `CMAKE_CXX_FLAGS` are **optional** and have the following effects:
- `-DSINGLE_PRECISION` selects single precision as the fundamental data type (default: `double`)
- `-DQRNG` ensures that the `QRNGRandomDevice`-type RNG is used. Provided without `-DQRNG_USE_DEVICE=1` this will result in a file-based RNG, which defaults to `/dev/random`.
- `-DQRNG_USE_DEVICE=1` *in conjunction with* `-DQRNG` selects the QRNG library functions as source of entropy, thereby electing to use a physical QRNG device. This option
without `-DQRNG` has no effect!

**Note**: Setting the `-DQRNG_USE_DEVICE` flag in `CMAKE_CXX_FLAGS` **will** override the setting made by CMake.
By default CMake will set `-DQRNG_USE_DEVICE=1` if the `lib` directory of the repository contains the QRNG library
(specifically `lib/QRNGlibrary/include/qrng_api.h`). Otherwise the falg will be set to 0 to ensure that
the code compiles when the QRNG library is not present. In such case the QRNG is replaced by `/dev/random` as entropy
source.

Currently there is no compile-time option to select the actual physical QRNG device used. This is accomplished by editing `QRNGRandomDevice.h` to select the 
appropriate initialisation parameters. The code was developedto use [Quantum Dice's](https://www.quantum-dice.com/) quantum random number generators.
It should work with any other QRNG/TRNG device by either modifying the `fillBuffer` function of the `QRNGRandomDevice` class or by using a file as a proxy.

## External libraries

## QD QRNG library

The Quantum Dice QRNG device library is optional and only required if one has access to their QRNG.
If you have access to such a device the library files should be placed in `lib/QRNGlibrary/lib`
and the headers in `lib/QRNGlibrary/include` for CMake to include them automatically.


## Tina's Random Number Generator library (trng)

In a parallel environment we elect to use a dedicated parallel pseudo-random number generator. 
Due to the [scalable parallel random number generator library](http://sprng.cs.fsu.edu/) (SPRNG) being outdated
we elected to use [Tina's random number generator library](https://www.numbercrunch.de/trng), unfortunately abbreviated
to TRNG. Specifically we use [v4.24 of TRNG](https://github.com/rabauke/trng4/tree/v4.24) and therein the parallel
PRNGs `mrg5`, `mrg5s` and `yarn5` following the rule-of-thumb often presented in introductions to MC simulations
of preferring PRNGs with long periods.

The source code expects the library to be located in the `lib/trng-library/lib64` directory and the 
library headers to be present in `lib/trng-library/include` directory.

**Note**: The functionality used in the present code has been tested with TRNG v4.15 and is likely to work
with older versions, too. Hence should the compilation of the newer version fail we suggest using an older version,
whose build requirements are less demanding.

### Compilation

To compile TRNG first initialise and update the submodule
`git submodule init`

`git submodule update`

Then create a `build` directory in `lib/trng` and step down into it:
`mkdir build && cd build`

Configure the library to compile without using CUDA and without additional tests. The latter is to prevent configuration failures
due to a missing [Catch2](https://github.com/catchorg/Catch2) test framework. Omission of the tests does not impact the functionality
of the code.
`cmake -DTRNG_ENABLE_TESTS=OFF`
`make && make install`


# Additional programs

Apart from the direct sampling Monte Carlo calculation of $\pi$ this repository also contains programs used
to evaluate the uniformity of the point distribution on $[0,1]^2$ as well as a program to dump random bits
to file. 

## Entropy storage

The directory `entropyToFile` contains source code for a small program that can be used to store
random bits to a binary file. The code utilises the `QRNGRandomDevice` class defined in the corresponding
header in `DSMCPI/include` and can hence be used to store either the random bits provided by Quantum Dice's
QRNG, the random bits provided by `/dev/random` (default behaviour) or any file that is used in the `qrng_init_proxy` function
in `QRNGRandomDevice.h`.

The intent of this code is to facilitate reproducibility with a true RNG by storing entropy to a file and
utilising said file for later simulations. The entropy file can be used with the same `QRNGRandomDevice` class
by substituting its absolute path for `/dev/random` within the aforementioned function.

## Unit (hyper)-cube sampling

The directory `nDCoordinates` contains the source for a program used to generate a point distribution on $[0,1]^d$.
This distribution can then be used to analyse the uniformity and homogeneity of the sampling of the $d$-dimensional
cube by the given RNG.

This code utilises the `QRNGRandomDevice` class to draw random bits either from a QRNG or, in absence of the
required libraries and headers, from `/dev/random`. If `-DQRNG` is not defined during the compilation process 
(c.f. compilation section above) the code will default to the Mersenne Twister PRNG `std::mt19937` introduced
with C++11.

The intent of this program is to generate points sampling the unit cube using various RNGs _in conjunction_
with the `uniform_real_distribution`. Tests performed on the points generated by this code are thus
evaluating the pairing of the RNG and the distribution! 
It is possible that a given RNG library provides methods to generate a random number in $[0,1)$ but it is not guaranteed to do so.
Additionally use of such methods would require validation of the conversion of random bits to the random number. To avoid
having to perform such a validation and to avoid introducing differences between PRNG and QRNG distributions that
would be solely due to the conversion method we elected to use the distribution functions of C++11.
This should ensure that any differences in the point distributions are solely due to the RNG.

# Evaluation of experiments

The directory `evaluation` and its sub-directories contain Jupyter notebooks used to process raw data obtained with
the Snakemake workflows and to produce the figures presented in the paper from processed data.

To run the notebooks following Python libraries should be installed:
- NumPy
- SciPy (+ stats)
- Matplotlib
- Numba
- tabulate
- Jupyter

The evaluations were performed with Python 3.8.2 as well as 3.11.
The paths to the data are provided relative to the location of each notebook. It is therefore advisable to launch Jupyter in the `evaluation` directory.

The analysis notebooks are provided in their evaluated form to ensure preservation of information even with 
the data missing.

**Note**:
1. Some sections of the `SciPaperVisualisation` notebook will require a sizeable amount of RAM 
and take a long time to evaluate. These sections are marked appropriately and can be omitted without
loss of information as they have not been included in the paper.
2. Not all figures present in these notebooks have been incorporated into the publication.