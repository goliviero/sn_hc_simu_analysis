#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(sn_hc_simu_analysis VERSION "0.0.1")

set(PROJECT_CXX_FLAGS "-Wall -Wextra -Wshadow -pedantic -Werror -Wfatal-errors")

# - Prepend flags to any preexisting ones, allowing user
#   to supply additional ones or override if required
set(CMAKE_CXX_FLAGS "${PROJECT_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

#-----------------------------------------------------------------------
# Falaise is the main external dependency, and we know it will additionally
# search for and provide compatible versions of:
# - Bayeux
# - Boost
# - GSL
# - ROOT
# - Geant4

set(Falaise_DIR /sps/nemo/scratch/golivier/software/Falaise/Falaise-install/lib64/cmake/Falaise-3.0.0/)
set(Bayeux_DIR /sps/nemo/scratch/golivier/software/Bayeux/Bayeux-install/lib64/cmake/Bayeux-3.0.0/)

find_package(Falaise 3.0.0 REQUIRED)


set(EXECUTABLE_OUTPUT_PATH BuildProducts/bin/)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#

file(GLOB_RECURSE source_files src/*.h)
file(GLOB_RECURSE source_files src/*.cxx)
file(GLOB_RECURSE source_files src/*.cpp)

add_executable(hc_sort_and_analysis_data source/hc_sort_and_analysis_data.cxx)

target_link_libraries(hc_sort_and_analysis_data Falaise::Falaise)
