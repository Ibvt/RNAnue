cmake_minimum_required (VERSION 3.4)
project (RNAnue VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

message("Configuring Config.h")
message("Input file: ${CMAKE_CURRENT_SOURCE_DIR}/include/Config.h.in")
message("Output file: ${CMAKE_CURRENT_BINARY_DIR}/include/Config.h")


# configure header file to pass the version number to  the source code
#configure_file(./include/Config.h.in ./include/Config.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/Config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/include/Config.h)

set(Boost_USE_MULTITHREADED TRUE)

set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/seqan3/build_system")
find_package (SeqAn3 REQUIRED)
find_package(OpenMP)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

set(Boost_DEBUG ON)
### boost
find_package(Boost 1.56.0 REQUIRED COMPONENTS
        program_options filesystem)
if(Boost_FOUND)
    include_directories( ${Boost_INCLUDE_DIRS} )
endif()

file(GLOB SOURCES "./source/*.cpp")
include_directories(./include)

add_executable (RNAnue ${SOURCES})
target_link_libraries (RNAnue seqan3::seqan3)
target_link_libraries (RNAnue Boost::program_options Boost::filesystem)
if(OpenMP_CXX_FOUND)
    target_link_libraries (RNAnue OpenMP::OpenMP_CXX)
endif()


# include RNAlib
#include(ExternalProject)
#ExternalProject_Add(project_RNAlib
#    URL https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz
#    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/RNAlib-2.4.14
#    CONFIGURE_COMMAND ""
#    BUILD_COMMAND make
#    INSTALL_COMMAND make install
#    PREFIX=${CMAKE_CURRENT_BINARY_DIR}/RNAlib-2-4.14
#)
#ExternalProject_Get_Property(project_RNAlib install_dir)
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#include_directories(${install_dir}/include/)
#add_dependencies(RNAnue project_RNAlib)
#target_link_libraries (RNAnue ${install_dir}/lib/RNAlib.a)

