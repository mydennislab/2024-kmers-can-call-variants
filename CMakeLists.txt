cmake_minimum_required (VERSION 3.4)
project (kDecontaminer C CXX)
set(default_build_type "Release")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC -lgomp -lrt -fopenmp -pedantic -O3 -Ofast")

add_subdirectory(lib/kmerDecoder)
add_subdirectory(lib/zstr)

include_directories("${PROJECT_SOURCE_DIR}/include")

include_directories("${PROJECT_SOURCE_DIR}/lib/kmerDecoder/lib/parallel-hashmap/")
include_directories("lib/progressbar/include")


CPMAddPackage(
        NAME kseq
        SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/lib/kmerDecoder/kseq
)
set(KSEQ_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/lib/kmerDecoder/lib/kseq/include")
include_directories($KSEQ_INCLUDE_DIRS)

include_directories("${PROJECT_SOURCE_DIR}/lib/argh")

# add_executable (peReadsStats reads_stats.cpp src/kseqReader.cpp)
# target_link_libraries (peReadsStats kProcessor z)
# target_include_directories(peReadsStats INTERFACE ${kProcessor_INCLUDE_PATH})

# add_executable (analyzeIdx analyze_index.cpp)
# target_link_libraries (analyzeIdx kProcessor)
# target_include_directories(analyzeIdx INTERFACE ${kProcessor_INCLUDE_PATH})

# add_executable (cDBG_partitioner cDBG_partitioner.cpp)
# target_link_libraries (cDBG_partitioner kProcessor z)
# target_include_directories(cDBG_partitioner INTERFACE ${kProcessor_INCLUDE_PATH})

# add_executable (commonKmers commonKmers.cpp)
# target_link_libraries (commonKmers kProcessor z)
# target_include_directories(commonKmers INTERFACE ${kProcessor_INCLUDE_PATH})

# add_executable (fastq_partitioner fastq_partitioner.cpp)
# target_link_libraries (fastq_partitioner kProcessor z)
# target_include_directories(fastq_partitioner INTERFACE ${kProcessor_INCLUDE_PATH})

#add_executable (testKseq test_kseq.cpp)
#target_link_libraries (testKseq z)



# genesIndexing executable (Have an issue after some time in the execution)
#add_executable (genesIndexing genes_indexing.cpp)
#target_link_libraries (genesIndexing kProcessor pthread z)
#target_include_directories(genesIndexing INTERFACE ${kProcessor_INCLUDE_PATH})


# 1) Convert the reference to bin
add_executable(refToBin 1_ref_to_bin.cpp)
target_link_libraries(refToBin kmerDecoder pthread z zstr::zstr)


add_executable(fqPart 2_fastq_partitioner.cpp src/helper_functions.cpp)
target_link_libraries(fqPart kmerDecoder pthread z zstr::zstr kseq::kseq)


# Testing
add_executable(test test/test_scoring.cpp src/helper_functions.cpp)
target_link_libraries(test)