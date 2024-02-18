// create inverted kmers index in parallel
#include "string"
#include <stdexcept>
#include "tuple"
#include <sys/stat.h>
#include <zlib.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cassert>
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"
#include <kmerDecoder.hpp>
#include <helper_functions.hpp>


using BINS_PHMAP = phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>,
    phmap::priv::hash_default_hash<std::string>,
    phmap::priv::hash_default_eq<std::string>,
    std::allocator<std::pair<const std::string, phmap::flat_hash_set<uint64_t>>>,
    1,
    std::mutex>;



int main(){
    
    kmerDecoder * KMER_HASHER = kmerDecoder::getInstance(KMERS, mumur_hasher, { {"kSize", 31} });
    KMER_HASHER->hash_kmer("ACTGCGTAGCTAGCTAGC");



    


}


