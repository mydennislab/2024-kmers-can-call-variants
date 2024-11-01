cmake_minimum_required(VERSION 3.4)
project(kDecontaminer C CXX)

# Set C++ standard and position-independent code
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
        set(default_build_type "Release")
        set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
endif()

# Common compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -O3 -Ofast")

# Find OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
        message(STATUS "Found OpenMP")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Add subdirectories
add_subdirectory(lib/kmerDecoder)
add_subdirectory(lib/zstr)

# Include directories
include_directories(
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/lib/kmerDecoder/lib/parallel-hashmap/"
        "${PROJECT_SOURCE_DIR}/lib/kmerDecoder/lib/kseq/include"
        "${PROJECT_SOURCE_DIR}/lib/argh"
)

# CPM Package for kseq
CPMAddPackage(
        NAME kseq
        SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/lib/kmerDecoder/kseq
)
set(KSEQ_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/lib/kmerDecoder/lib/kseq/include")
include_directories(${KSEQ_INCLUDE_DIRS})

# Add executable
add_executable(refToBin 1_ref_to_bin.cpp)

# Link libraries
target_link_libraries(refToBin
        kmerDecoder
        pthread
        z
        zstr::zstr
        OpenMP::OpenMP_CXX
)

# Link 'rt' library only on Linux
if(UNIX AND NOT APPLE)
        target_link_libraries(refToBin rt)
endif()
