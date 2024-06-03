# CMake version check
cmake_minimum_required(VERSION 3.13)  

# Create project
project(dalitz_studies)

# Enable c++17 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the source files
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Include the header files directory
include_directories("framework")

# Add eigen
include_directories("$ENV{CONDA_PREFIX}/include/eigen3")

# Add qft++
include_directories("$ENV{CONDA_PREFIX}/external/qft/include")
link_directories("$ENV{CONDA_PREFIX}/external/qft/lib")

# Add Boost
find_package(Boost REQUIRED COMPONENTS program_options)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION = ${Boost_VERSION}")
else()
    message(FATAL_ERROR "Could not find Boost")
endif()

# Add ROOT
find_package(ROOT REQUIRED COMPONENTS Core Imt RIO Net Hist Graf Graf3d Gpad ROOTVecOps Tree TreePlayer Rint Postscript Matrix Physics MathCore Thread MultiProc ROOTDataFrame Minuit2 MathMore)
if(NOT ROOT_FOUND)
    message(FATAL_ERROR "Could not find ROOT")
endif()


include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

# Create an executable for each source file
foreach(SOURCE_FILE ${SOURCES})
    # Get the file name without directory and extension
    get_filename_component(EXECUTABLE_NAME ${SOURCE_FILE} NAME_WE)

    # Add the executable
    add_executable(${EXECUTABLE_NAME} ${SOURCE_FILE})

    # Add the CONDA and ROOT compiler flags
    # target_compile_options(${EXECUTABLE_NAME} PRIVATE ${ROOT_CFLAGS})

    # Link the library from CONDA_PREFIX with the appropriate extension
    target_link_libraries(${EXECUTABLE_NAME} qft++ ${Boost_LIBRARIES} ${ROOT_LIBRARIES})

    # Set the output directory for the binaries
    set_target_properties(${EXECUTABLE_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
    )
endforeach()

# If you want the binaries in ./bin instead of the default CMake binary directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
