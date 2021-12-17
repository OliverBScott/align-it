cmake_minimum_required(VERSION 3.14)
project(align_it)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

option(BUILD_RDKIT_SUPPORT "Build against the RDKit instead of OpenBabel3" OFF)
option(BUILD_PYTHON_SUPPORT "Build Python wrappers (requires RDKit support)" OFF)

if(BUILD_PYTHON_SUPPORT AND (NOT BUILD_RDKIT_SUPPORT))
    message(WARNING "Disabling BUILD_PYTHON_SUPPORT because BUILD_RDKIT_SUPPORT is not enabled")
    set(BUILD_PYTHON_SUPPORT OFF)
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
set(ALIGNIT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 4)
file(WRITE  ${ALIGNIT_INCLUDE_DIR}/config.h "#define ALIGNIT_VERSION ${VERSION_MAJOR}\n")
file(APPEND ${ALIGNIT_INCLUDE_DIR}/config.h "#define ALIGNIT_RELEASE ${VERSION_MINOR}\n")
file(APPEND ${ALIGNIT_INCLUDE_DIR}/config.h "#define ALIGNIT_SUBRELEASE ${VERSION_PATCH}\n")

set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0)
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)

include_directories(${ALIGNIT_INCLUDE_DIR})
if(NOT BUILD_RDKIT_SUPPORT)
    message(STATUS "Building Align-it against OpenBabel3")
    find_package(OpenBabel3 REQUIRED)
    include_directories(${OPENBABEL3_INCLUDE_DIR})
else(NOT BUILD_RDKIT_SUPPORT)
    message(STATUS "Building Align-it against RDKit")
    add_definitions(-DUSE_RDKIT)
    find_package(RDKit REQUIRED)
    find_package(Threads)
    set(RDKit_THREAD_LIBS Threads::Threads)
    find_package(Boost REQUIRED)
endif(NOT BUILD_RDKIT_SUPPORT)

add_library(alignit_lib SHARED src/siMath.cpp src/coordinate.cpp src/pharmacophore.cpp
    src/stringTokenizer.cpp src/chargeFuncCalc.cpp src/solutionInfo.cpp src/utilities.cpp
    src/pharMerger.cpp src/hybridCalc.cpp src/alignment.cpp src/functionMapping.cpp
    src/result.cpp src/addBest.cpp src/compScore.cpp src/addBest.cpp
    src/aromFuncCalc.cpp src/hAccFuncCalc.cpp src/hDonFuncCalc.cpp src/lipoFuncCalc.cpp
    src/calcPharm.cpp src/alignLib.cpp src/mainWar.cpp
)
target_include_directories(alignit_lib PUBLIC ${OPENBABEL3_INCLUDE_DIR} ${RDKIT_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
target_link_libraries(alignit_lib ${OPENBABEL3_LIBRARIES} ${RDKIT_LIBRARIES} ${RDKit_THREAD_LIBS})
install(TARGETS alignit_lib DESTINATION lib)

add_executable(${PROJECT_NAME} src/options.cpp src/main.cpp src/mainErr.cpp src/printHeader.cpp
    src/printInfo.cpp src/printUsage.cpp src/printProgress.cpp src/parseCommandLine.cpp
    src/logOut.cpp src/logScores.cpp src/logPharmacophores.cpp src/getExt.cpp
)
target_link_libraries(${PROJECT_NAME} alignit_lib ${OPENBABEL3_LIBRARIES} ${RDKIT_LIBRARIES} ${RDKit_THREAD_LIBS})
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
message("Installation directory of executable: ${CMAKE_INSTALL_PREFIX}/bin")

if(BUILD_PYTHON_SUPPORT)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_subdirectory(src/Wrap)
    message(STATUS "Installing to: ${Python3_SITEARCH}/pyalignit")
    install(DIRECTORY pyalignit DESTINATION ${Python3_SITEARCH})
endif(BUILD_PYTHON_SUPPORT)
