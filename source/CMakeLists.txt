cmake_minimum_required (VERSION 3.4)
project (RNAnue VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# configure header file to pass the version number to the source code
configure_file(../include/Config.h.in ../include/Config.h)

set(Boost_USE_MULTITHREADED TRUE)

set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../../seqan3/build_system")
find_package (SeqAn3 REQUIRED)

#find_package(OpenMP)

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

file(GLOB SOURCES "*.cpp")
include_directories(../include)

add_executable (RNAnue ${SOURCES} /usr/local/include/ViennaRNA/)
target_link_libraries (RNAnue seqan3::seqan3)
target_link_libraries (RNAnue Boost::program_options Boost::filesystem)

#target_include_directories(RNAnue PUBLIC "${PROJECT_BINARY_DIR}")
#target_link_libraries (RNAnue OpenMP::OpenMP_CXX)
