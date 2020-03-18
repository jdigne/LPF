cmake_minimum_required(VERSION 2.6)
project(lpf)

find_package(OpenMP)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
endif()
if(UNIX)
    SET(CMAKE_CXX_FLAGS_RELEASE
        "-O3 -Wall -Wextra -std=c++11 -Wfatal-errors ${OpenMP_CXX_FLAGS}")
    SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g  -Wall -Wextra -std=c++11 -Wfatal-errors")
else()
    set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
endif()


#message(${CMAKE_CXX_COMPILER})
#message(${CMAKE_CXX_FLAGS_RELEASE})
#message(${CMAKE_CXX_FLAGS_DEBUG})

add_executable(lpf src/main.cpp
    src/LPF.cpp
    src/io.cpp
    src/Parameterizer.cpp
    src/DescriptorComputation.cpp
    src/DescriptorAnalysis.cpp
    src/DescriptorReconstruction.cpp
    src/Denoising.cpp
    src/Optimization.cpp)

install(TARGETS lpf RUNTIME DESTINATION bin)



#eigen directory
include_directories(./eigen3/)
