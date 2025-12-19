cmake_minimum_required(VERSION 3.18)

project(smartsed LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Options
# ------------------------------------------------------------------------------

option(ENABLE_CUDA "Enable CUDA support" OFF)

# ------------------------------------------------------------------------------
# Build type (Debug / Release)
# ------------------------------------------------------------------------------

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build: Debug or Release" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CUDA enabled: ${ENABLE_CUDA}")

# ------------------------------------------------------------------------------
# C++ standard
# ------------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------------------------------------------------------------------
# MPI
# ------------------------------------------------------------------------------

find_package(MPI REQUIRED)

# ------------------------------------------------------------------------------
# Eigen in case no CUDA
# ------------------------------------------------------------------------------

if(NOT ENABLE_CUDA)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    message(STATUS "Eigen found: ${Eigen3_VERSION}")
endif()

# ------------------------------------------------------------------------------
# CUDA (optional)
# ------------------------------------------------------------------------------

if(ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)

    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

# ------------------------------------------------------------------------------
# Executable
# ------------------------------------------------------------------------------

add_executable(main.exe
    main_final_H.cpp
    utils_H.cpp
    cuda_utils_loop_H.cu
)

# ------------------------------------------------------------------------------
# Include directories
# ------------------------------------------------------------------------------

target_include_directories(main.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/eigen_forma
    ${CMAKE_SOURCE_DIR}/include/eigen_forma/Utilities
)

# ------------------------------------------------------------------------------
# Compile definitions
# ------------------------------------------------------------------------------

target_compile_definitions(main.exe PRIVATE
    ENABLE_TIMING
)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_definitions(main.exe PRIVATE NDEBUG)
endif()

if(ENABLE_CUDA)
    target_compile_definitions(main.exe PRIVATE ENABLE_CUDA)
endif()

# ------------------------------------------------------------------------------
# Compile options
# ------------------------------------------------------------------------------

target_compile_options(main.exe PRIVATE
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>:-O3>
)

# ------------------------------------------------------------------------------
# Link libraries
# ------------------------------------------------------------------------------

target_link_libraries(main.exe PRIVATE
    MPI::MPI_CXX
)

if(NOT ENABLE_CUDA)
    target_link_libraries(main.exe PRIVATE Eigen3::Eigen)
endif()


if(ENABLE_CUDA)
    target_link_libraries(main.exe PRIVATE CUDA::cudart)
endif()

