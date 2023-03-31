#include <iostream>
#include <cstdint>
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include <ctime>
#include <omp.h>
#include <glob.h>
#include <string>
#include <stdexcept>
#include "parallel_hashmap/phmap_dump.h"
#include <cstdlib>
#include "zstr.hpp"
#include <kmerDecoder.hpp>

using namespace std;
// using namespace phmap;

typedef std::chrono::high_resolution_clock Time;

string reformat_long_number(uint64_t number) {
    string number_str = to_string(number);
    string new_number_str = "";
    int counter = 0;
    for (int i = number_str.length() - 1; i >= 0; i--) {
        new_number_str = number_str[i] + new_number_str;
        counter++;
        if (counter == 3 && i != 0) {
            new_number_str = "," + new_number_str;
            counter = 0;
        }
    }
    return new_number_str;
}

int main(int argc, char** argv) {

    if (argc != 4) {
        cout << "run: ./refToBin <fasta_path> <kSize> <output_path>" << endl;
        exit(1);
    }


    string fasta_path = argv[1];
    int kSize = stoi(argv[2]);
    string output_path = argv[3];
    int chunk_size = 1000;


    auto begin_time = Time::now();

    // phmap::flat_hash_set<uint64_t> hashes;
    phmap::flat_hash_map<uint64_t, uint32_t> hashes;


    std::string base_filename = fasta_path.substr(fasta_path.find_last_of("/\\") + 1);
    base_filename = base_filename.substr(0, base_filename.find('_'));

    kmerDecoder* REF_KMERS = kmerDecoder::getInstance(fasta_path, chunk_size, KMERS, mumur_hasher, { {"kSize", kSize} });
    uint64_t total_kmers = 0;

    while (!REF_KMERS->end() && !REF_KMERS->end()) {
        REF_KMERS->next_chunk();

        flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1 = REF_KMERS->getKmers()->begin();
        flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1_end = REF_KMERS->getKmers()->end();

        while (seq1 != seq1_end) {
            for (auto const kRow : seq1->second) {
                hashes[kRow.hash]++;
                total_kmers++;
            }
            seq1++;
        }
    }

    cout << "inserted " << reformat_long_number(hashes.size()) << " unique kmers out of " << reformat_long_number(total_kmers) << " ." << endl;
    string out_path = output_path;
    phmap::BinaryOutputArchive ar_out(out_path.c_str());
    hashes.phmap_dump(ar_out);
    cout << "Conversion done in " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;
}