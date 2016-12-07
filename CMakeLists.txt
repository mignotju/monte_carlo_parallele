cmake_minimum_required(VERSION 2.8)
project(mc-pricer C CXX)

# Prefix
if (NOT PREFIX)
    set(PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif (NOT PREFIX)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE STRING "Installation prefix." FORCE)

if (CXX)
    set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()

# Detect MPI
set(MPI_LIBS)
find_package(MPI REQUIRED)
if (MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set (MPI_LIBS ${MPI_LIBS} ${MPI_C_LIBRARIES})
    add_definitions(-DOMPI_SKIP_MPICXX)
else (MPI_CXX_FOUND)
    message(STATUS "MPI not found.")
endif (MPI_CXX_FOUND)

# Detect OMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS ${CMAKE_C_CFLAGS} ${OpenMP_C_FLAGS})
    set (CMAKE_CXX_FLAGS ${CMAKE_CXX_CFLAGS} ${OpenMP_CXX_FLAGS})
else (OPENMP_FOUND)
    message(STATUS "OpenMP not found.")
endif (OPENMP_FOUND)

# Detect Boost
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
message(STATUS "Boost include directory ${Boost_INCLUDE_DIR}")

# Detect PNL
find_package(Pnl REQUIRED)
set(LIBS ${LIBS} ${PNL_LIBRARIES})
include_directories(${PNL_INCLUDE_DIRS})
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
endif()

# Add math lib
set(LIBS ${LIBS} m)

# Testing the parser
#add_executable(parser-test src/parser.cpp src/test_parser.cpp)
#target_link_libraries(parser-test ${LIBS})
#add_executable(mc-example src/mc_example.cpp)
#target_link_libraries(mc-example ${LIBS})



#tests des main
#add_executable(main0 src/MonteCarlo.cpp src/OptionBasket.cpp src/OptionAsiatique.cpp src/OptionPerformance.cpp src/parser.cpp src/BlackScholesModel.cpp src/main0.cpp)
#target_link_libraries(main0 ${LIBS})
#add_executable(main_t src/MonteCarlo.cpp src/OptionBasket.cpp src/parser.cpp src/BlackScholesModel.cpp src/main_t.cpp)
#target_link_libraries(main_t ${LIBS})



# Test Payoff Option Basket
#add_executable(test_basket tests/test_basket.cpp src/OptionBasket.cpp)
#target_link_libraries(test_basket ${LIBS})

# Test Payoff Option Asiatique 
#add_executable(test_asiatique tests/test_asiatique.cpp src/OptionAsiatique.cpp)
#target_link_libraries(test_asiatique ${LIBS})

# Test Payoff Option Performance 
#add_executable(test_performance tests/test_performance.cpp src/OptionPerformance.cpp)
#target_link_libraries(test_performance ${LIBS})

# Test MonteCarlo 
#add_executable(test_MonteCarlo src/parser.cpp tests/test_MonteCarlo.cpp src/MonteCarlo.cpp src/OptionBasket.cpp src/OptionPerformance.cpp src/OptionAsiatique.cpp src/BlackScholesModel.cpp)
#target_link_libraries(test_MonteCarlo ${LIBS})

# Test Simulation 
#add_executable(test_asset src/parser.cpp tests/test_asset.cpp src/BlackScholesModel.cpp)
#target_link_libraries(test_asset ${LIBS})

# Test Delta
#add_executable(test_delta2 tests/test_delta2.cpp src/parser.cpp src/BlackScholesModel.cpp src/MonteCarlo.cpp src/OptionBasket.cpp src/OptionPerformance.cpp src/OptionAsiatique.cpp)
#target_link_libraries(test_delta2 ${LIBS})

# add_executable(mc-pricer src/parser.cpp src/pricer.cpp)
# target_link_libraries(mc-pricer ${LIBS})
# install(TARGETS mc-pricer DESTINATION bin)

# Test Simulation 
#add_executable(test_delta tests/test_delta.cpp src/MonteCarlo.cpp src/parser.cpp src/BlackScholesModel.cpp src/OptionBasket.cpp src/OptionAsiatique.cpp src/OptionPerformance.cpp)
#target_link_libraries(test_delta ${LIBS})
#add_executable(test_couverture tests/test_couverture.cpp src/parser.cpp src/MonteCarlo.cpp src/Simulation.cpp src/BlackScholesModel.cpp src/OptionBasket.cpp src/OptionAsiatique.cpp src/OptionPerformance.cpp)
#target_link_libraries(test_couverture ${LIBS})
#add_executable(test_simul_market tests/test_simulMarket.cpp src/parser.cpp src/MonteCarlo.cpp src/Simulation.cpp src/BlackScholesModel.cpp src/OptionBasket.cpp src/OptionAsiatique.cpp src/OptionPerformance.cpp)
#target_link_libraries(test_simul_market ${LIBS})

# Executable final : pricer
add_executable(pricer src/pricer.cpp src/parser.cpp src/MonteCarlo.cpp src/Simulation.cpp src/BlackScholesModel.cpp src/OptionBasket.cpp src/OptionAsiatique.cpp src/OptionPerformance.cpp)
target_link_libraries(pricer ${LIBS})

# Print compiler flags
get_directory_property(DirDefs COMPILE_DEFINITIONS)
message(STATUS "COMPILE_DEFINITIONS = ${DirDefs}")
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}" )
    message(STATUS "C Flags: ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}" )
else()
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}" )
    message(STATUS "C Flags: ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}" )
endif()
