
Exploration of effects of entropy (random numbers) sources on Monte Carlo simulations.

# Direct Sampling Monte Carlo estimation of Pi (DSMCPI)

The `pi-estimation` branch of this repository contains the source
code, [Snakemake](https://snakemake.readthedocs.io/en/stable) workflows and Jupyter notebooks
used to test the effects of different random number generators (RNGs) on simple
Monte Carlo applications.

This branch also contains
a tool to dump random bits from the QRNG device (or `/dev/random`) to a binary file and
a tool to convert those binary files into lists of n-dimensional vectors in $[0,1]^d$.

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
