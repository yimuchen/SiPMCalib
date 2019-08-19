cmake_minimum_required(VERSION 3.0)
project(SiPMCalibAnalyze)

## Compilter settings
set(CMAKE_CXX_COMPILE_FEATURES cxx_range_for )
set(CMAKE_CXX_COMPILE_FEATURES cxx_variadic_templates)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options("-O3")
add_compile_options("-Wall")
add_compile_options("-Wno-undef")

# Command package setup
find_package(PythonLibs   3 REQUIRED)
find_package(PythonInterp 3 REQUIRED)
find_package( Boost REQUIRED
  COMPONENTS program_options python3 )
find_package( ROOT REQUIRED
  COMPONENTS RooFitCore RooFit Spectrum)
find_package( GSL REQUIRED )


# Declaring make directory
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/lib )
include_directories( ${CMAKE_HOME_DIRECTORY}/include )
include_directories( ${ROOT_INCLUDE_DIRS} )
include_directories( /home/ensc/Homework/CMSSW/ )
include_directories( /home/ensc/Projects/UserUtils/include/ )
link_directories( /home/ensc/Projects/UserUtils/lib/ )

## Libraries compile instructions are defined in directory
add_subdirectory(SiPMCalc)
