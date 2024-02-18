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
#include <queue>
#include <kSpider.hpp>


using BINS_MAP = phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>,
    phmap::priv::hash_default_hash<std::string>,
    phmap::priv::hash_default_eq<std::string>,
    std::allocator<std::pair<const std::string, phmap::flat_hash_set<uint64_t>>>,
    12,
    std::mutex
>;
using LEGENDS_MAP = phmap::parallel_flat_hash_map<uint64_t,
    std::vector<uint32_t>,
    std::hash<uint64_t>,
    std::equal_to<uint64_t>,
    std::allocator<std::pair<const uint64_t, vector<uint32_t>>>,
    4>; // 6 submaps because colors will grow

using LEGENDS_MAP_OLD = phmap::parallel_flat_hash_map<uint64_t, std::vector<uint32_t>>;

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

// thanks to https://stackoverflow.com/a/8615450/3371177
inline std::vector<std::string> glob2(const std::string& pattern) {
    using namespace std;

    // glob struct resides on the stack
    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    // do the glob operation
    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if (return_value != 0) {
        globfree(&glob_result);
        stringstream ss;
        ss << "glob() failed with return_value " << return_value << endl;
        throw std::runtime_error(ss.str());
    }

    // collect all the filenames into a std::list<std::string>
    vector<string> filenames;
    for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
        filenames.push_back(string(glob_result.gl_pathv[i]));
    }

    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}


int main(int argc, char** argv) {

    if (argc != 4) {
        cout << "run: ./refToBin <kSize> <min_abundance> <output_file>" << endl;
        exit(1);
    }


    string fasta_path = argv[1];
    int kSize = stoi(argv[2]);
    string output_path = argv[3];
    int chunk_size = 1000;

    kSpider::bins_indexing("", kSize, "km", 1000000, 222);
}