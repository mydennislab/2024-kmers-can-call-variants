// create inverted kmers index in parallel
#include "string"
#include "kseq/kseq.h"
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

// KSEQ_INIT(gzFile, gzread)

using namespace std;
using namespace phmap;


class fastqWriter {

public:
    ofstream fileStream_r1, fileStream_r2;

    explicit fastqWriter(string& filename_prefix) {
        this->fileStream_r1.open(filename_prefix + "_R1.fastq");
        this->fileStream_r2.open(filename_prefix + "_R2.fastq");
    }

    void write(kseq_t* seq1, kseq_t* seq2) {

        this->fileStream_r1 << '@';
        this->fileStream_r1 << seq1->name.s;
        this->fileStream_r1 << ' ';
        this->fileStream_r1 << seq1->comment.s;
        this->fileStream_r1 << endl;
        this->fileStream_r1 << seq1->seq.s;
        this->fileStream_r1 << "\n+\n";
        this->fileStream_r1 << seq1->qual.s;
        this->fileStream_r1 << endl;

        this->fileStream_r2 << '@';
        this->fileStream_r2 << seq2->name.s;
        this->fileStream_r2 << ' ';
        this->fileStream_r2 << seq2->comment.s;
        this->fileStream_r2 << endl;
        this->fileStream_r2 << seq2->seq.s;
        this->fileStream_r2 << "\n+\n";
        this->fileStream_r2 << seq2->qual.s;
        this->fileStream_r2 << endl;
    }

    void close() {
        fileStream_r1.close();
        fileStream_r2.close();
    }

};

inline bool file_exists(const std::string& name) {
    struct stat buffer
    {
    };
    return (stat(name.c_str(), &buffer) == 0);
}

template <typename KeyType, typename ValueType>
int count_intersections(const std::vector<KeyType>& s, const flat_hash_map<KeyType, ValueType>& m) {
    int intersection_count = 0;

    for (const auto& key : s) {
        if (m.find(key) != m.end()) {
            ++intersection_count;
        }
    }

    return intersection_count;
}


int main() {

    // KMER_HASHER->hash_kmer("ACTGCGTAGCTAGCTAGC");


    /*
        1- Load all bins
    */

    BINS_PHMAP bins;
    load_all_bins("/home/mabuelanin/dib-dev/2023-decontamination/data/bins", &bins, 3);
    cout << "All genomes has been loaded" << endl;

    // Mapping genomes to ids
    flat_hash_map<string, int> genome_to_id;
    flat_hash_map<int, string> id_to_genome;
    int id_counter = 0;
    for (const auto& [genome_name, _kmers] : bins) {
        genome_to_id[genome_name] = ++id_counter;
        id_to_genome[id_counter] = genome_name;
    }

    /*
        2- Creating file handlers
    */

    string PE_1_reads_file = "/home/mabuelanin/dib-dev/2023-decontamination/data/samples/subset_Ast25B_R1_001.fastq.gz";
    string PE_2_reads_file = "/home/mabuelanin/dib-dev/2023-decontamination/data/samples/subset_Ast25B_R2_001.fastq.gz";
    int chunk_size = 10000;
    int kSize = 21;
    kmerDecoder* KMER_HASHER = kmerDecoder::getInstance(KMERS, mumur_hasher, { {"kSize", kSize} });




    if (!file_exists(PE_1_reads_file)) {
        throw std::runtime_error("Could not open the unitigs fasta file");
    }

    if (!file_exists(PE_2_reads_file)) {
        throw std::runtime_error("Could not open the unitigs fasta file");
    }


    std::string base_filename = PE_1_reads_file.substr(PE_1_reads_file.find_last_of("/\\") + 1);
    base_filename = base_filename.substr(0, base_filename.find('R1'));

    map<string, fastqWriter*> fasta_writer;

    cerr << "Creating fasta file handlers" << endl;
    for (auto& [genome_name, _] : bins) {
        string file_name = "genome_" + genome_name + "_readsPartition";
        fasta_writer[genome_name] = new fastqWriter(file_name);
    }

    string unmapped_file_name = base_filename + "_unmapped_partition.fa";
    fasta_writer["unmapped"] = new fastqWriter(unmapped_file_name);


    /*
        3- Processing reads
    */

    gzFile fp_1, fp_2;
    kseq_t* kseq_1, * kseq_2;
    fp_1 = gzopen(PE_1_reads_file.c_str(), "r");
    fp_2 = gzopen(PE_2_reads_file.c_str(), "r");

    kseq_1 = kseq_init(fp_1);
    kseq_2 = kseq_init(fp_2);

    int Reads_chunks_counter = 0;

    cout << "Partitioning the reads" << endl;
    for (int seqCounter = 0; kseq_read(kseq_1) >= 0 && kseq_read(kseq_2) >= 0; seqCounter++) {

        uint32_t seq_1_length = string(kseq_1->seq.s).size();
        uint32_t seq_2_length = string(kseq_2->seq.s).size();

        if (seq_1_length < kSize || seq_2_length < kSize) continue;

        std::string seq = string(kseq_1->seq.s) + string(kseq_2->seq.s);

        // Convert sequence to kmers
        vector<kmer_row> kmers;
        KMER_HASHER->seq_to_kmers(seq, kmers);
        // Convert kmer_row kmers to uint64_t hashes.
        vector<uint64_t> kmers_vec;
        kmers_vec.reserve(kmers.size()); // Preallocate memory        
        for (const auto& kmer : kmers) {
            kmers_vec.emplace_back(kmer.hash);
        }





        // iterate over kmers and get the genomes that contain them


        vector<uint32_t> kmers_matches;
        for (auto& [genome_name, genome_kmers] : bins) {
            
            size_t intersection_kmers = count_intersections(kmers_vec, genome_kmers);
            
            for (size_t i = 0; i < intersection_kmers; i++) {
                kmers_matches.emplace_back(genome_to_id[genome_name]);
            }
        }


        auto category = classify_and_match_read_kmers(kmers_matches, 0.1, 2.0);
        if (category.first == "unmapped") {
            fasta_writer["unmapped"]->write(kseq_1, kseq_2);
        }
        else if (category.first == "unique") {
            assert(category.second.size() == 1);
            fasta_writer[id_to_genome[category.second[0]]]->write(kseq_1, kseq_2);
        }
        else if (category.first == "ambiguous") {
            for (auto const& genomeID : category.second) {
                fasta_writer[id_to_genome[genomeID]]->write(kseq_1, kseq_2);
            }
        }


    }

}



// kmerDecoder * KMER_HASHER = kmerDecoder::getInstance(KMERS, mumur_hasher, { {"kSize", 31} });
//     KMER_HASHER->hash_kmer("ACTGCGTAGCTAGCTAGC");

