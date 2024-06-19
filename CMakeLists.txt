cmake_minimum_required(VERSION 3.22.1)
project(RNAnue VERSION 0.2.3)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS -fopenmp)

include(CMakePrintHelpers)

# configure header file to pass the version number to the source code
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/Config.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/include/Config.hpp)

###### SeqAn #####
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/seqan3/build_system")
find_package (seqan3 3.0 REQUIRED)
find_package(OpenMP)

###### Boost ######
set(Boost_DEBUG ON)
find_package(Boost 1.56.0 REQUIRED COMPONENTS
        program_options filesystem unit_test_framework)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

###### HTSlib ######
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED IMPORTED_TARGET htslib)
include_directories(${HTSLIB_INCLUDE_DIRS})

file(GLOB SOURCES "src/*.cpp")
add_executable(RNAnue ${SOURCES})
target_include_directories(RNAnue PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(RNAnue seqan3::seqan3)
target_link_libraries(RNAnue Boost::program_options Boost::filesystem)
target_link_libraries(RNAnue PkgConfig::HTSLIB)
#target_link_libraries(RNAnue -L/usr/local/include/ -lRNA)

cmake_print_properties(TARGETS RNAnue PROPERTIES TARGET_INCLUDE_DIRECTORIES)

###### Tests ######
file(GLOB TEST_SOURCES "test/*.cpp")
file(GLOB SOURCES "src/*.cpp")
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

add_executable(RNAnue_tests ${TEST_SOURCES} ${SOURCES})
target_include_directories(RNAnue_tests PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(RNAnue_tests seqan3::seqan3)
target_link_libraries(RNAnue_tests Boost::program_options Boost::filesystem Boost::unit_test_framework)
target_link_libraries(RNAnue_tests PkgConfig::HTSLIB)

