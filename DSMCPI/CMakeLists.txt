# CMAKE File for the DSMC Pi calculation with various RNGs
cmake_minimum_required(VERSION 3.20)

# Declare project but don't process compilers yet
project(DSMCPI NONE)

# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_CXX_STANDARD_REQUIRED C++11)
# End of setup and error checking.

# Now enable the compiler language support
enable_language(CXX)
enable_language(C)

#add BOOST dependency, requires CXX enabled
find_package(Boost COMPONENTS program_options REQUIRED)
#add OpenMP as requirement
find_package(OpenMP)

#Check if the QRNG API header is present
# if not we disable inclusion of QRNG sources
#and ensure that only proxy functions are used
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../lib/QRNGlibrary/include/qrng_api.h)
    add_compile_definitions(-DQRNG_USE_DEVICE=1)
    # Adding PCAP dependency for the QRNG library
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/")
    include(FindPCAP.cmake)
    find_package(PCAP REQUIRED)
    
    #add the sources for the QRNG library 
    include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/../lib/QRNGlibrary/include)
    #set preprocessor flag to compile QRNG functions
else()
    add_compile_definitions(-DQRNG_USE_DEVICE=0)
endif()


#add header directories
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#TRNG sources
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/../lib/trng-library/include)
#Collect source files from the source directory and store them in the `sources` variable
file(GLOB sources "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

#Define the name of the executable
add_executable(DSMCPI ${sources})

#link libraries
target_link_libraries(DSMCPI PUBLIC Boost::program_options)

if(OpenMP_CXX_FOUND)
        target_link_libraries(DSMCPI PUBLIC OpenMP::OpenMP_CXX)
else()
  message( FATAL_ERROR "Your compiler does not support OpenMP" )
endif()
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../lib/QRNGlibrary/include/qrng_api.h)
    message( INFO "Linking to QRNG library")
    target_link_libraries(DSMCPI PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../lib/QRNGlibrary/lib/libqrnglib.a)
    message( INFO "Linking to " ${PCAP_LIBRARY})
    target_link_libraries(DSMCPI PUBLIC ${PCAP_LIBRARY})
endif()
#add TRNG
message( INFO "Linking to parallel PRNG library: Tina's Random Number Generator Library")
target_link_libraries(DSMCPI PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../lib/trng-library/lib64/libtrng4.a)