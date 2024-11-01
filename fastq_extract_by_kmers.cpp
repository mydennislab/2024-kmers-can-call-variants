// create inverted kmers index in parallel
#include <string>
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
#include <argh.h>

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

void print_help()
{
    std::cout << "Usage: program_name [options]\n"
        "Options:\n"
        "  --R1 <file>          R1 input file (required)\n"
        "  --R2 <file>          R2 input file (required)\n"
        "  --genomes_dir <dir>  dir containing genomes bins created by refToBin\n"
        "  --kSize <int>    K-mer size (required)\n"
        "  --threads <number>   Number of genomes loading threads (default: 1)\n"
        "  --help, -h           Show this help message and exit\n";
}




int main(int argc, char** argv) {

    argh::parser cmdl;
    cmdl.add_params({ "--genomes_dir", "--threads", "--R1", "--R2", "--kSize", "--output" });
    cmdl.parse(argc, argv);

    if (cmdl[{"--help", "-h"}])
    {
        print_help();
        return 0;
    }

    int threads = 1; // Default value
    std::string R1_file;
    std::string R2_file;
    std::string output = "NA"; // Default value
    int kSize;
    std::string genomes_dir;

    // Check required arguments and retrieve their values
    if (cmdl({ "--threads" }, 1)) cmdl({ "--threads" }, 1) >> threads;
    
    if (cmdl({ "--R1" }))
        cmdl({ "--R1" }) >> R1_file;
    else
    {
        std::cerr << "Error: R1 file argument is required.\n";
        return 1;
    }

    if (cmdl({ "--R2" }))
        cmdl({ "--R2" }) >> R2_file;
    else
    {
        std::cerr << "Error: R2 file argument is required.\n";
        return 1;
    }

    if (cmdl({ "--genomes_dir" }))
        cmdl({ "--genomes_dir" }) >> genomes_dir;
    else
    {
        std::cerr << "Error: genomes_dir argument is required.\n";
        return 1;
    }

    if (cmdl({ "--kSize" }))
        cmdl({ "--kSize" }) >> kSize;
    else
    {
        std::cerr << "Error: kSize argument is required.\n";
        return 1;
    }

    if (cmdl({ "--output" }))
        cmdl({ "--output" }) >> output;

    std::cout << "Threads: " << threads << std::endl;
    std::cout << "R1 File: " << R1_file << std::endl;
    std::cout << "R2 File: " << R2_file << std::endl;
    std::cout << "Output: " << output << std::endl;
    std::cout << "K-mer size: " << kSize << std::endl;
    std::cout << "Genomes dir: " << genomes_dir << std::endl;
    std::cout << "----------------------------------------\n" << std::endl;
    // print line
    std::cout << "----------------------------------------\n" << std::endl;


    /*
        0- stats
    */

    Stats stats;
    flat_hash_map<string, Stats> sample_stats;


    /*
        1- Load all bins
    */

    BINS_PHMAP bins;
    load_all_bins(genomes_dir, &bins, threads);
    cout << "All genomes has been loaded" << endl;


    // Mapping genomes to ids
    flat_hash_map<string, int> genome_to_id;
    flat_hash_map<int, string> id_to_genome;
    int id_counter = 0;
    for (const auto& [genome_name, _kmers] : bins) {
        genome_to_id[genome_name] = ++id_counter;
        id_to_genome[id_counter] = genome_name;
        cout << genome_name << " -> " << id_counter << endl;
        // length
        cout << "Length: " << _kmers.size() << endl;
    }
    stats.set_id_to_genome(id_to_genome);

    /*
        2- Creating file handlers
    */

    kmerDecoder* KMER_HASHER = kmerDecoder::getInstance(KMERS, mumur_hasher, { {"kSize", kSize} });
    bool write_fastq = write;


    if (!file_exists(R1_file)) {
        throw std::runtime_error("Could not open the unitigs fasta file");
    }

    if (!file_exists(R2_file)) {
        throw std::runtime_error("Could not open the unitigs fasta file");
    }


    std::string base_filename = R1_file.substr(R1_file.find_last_of("/\\") + 1);
    base_filename = base_filename.substr(0, base_filename.find('R1'));

    map<string, fastqWriter*> fasta_writer;
    if (write_fastq) {
        cerr << "Creating fasta file handlers" << endl;
        for (auto& [genome_name, _] : bins) {
            string file_name = "genome_" + genome_name + "_readsPartition";
            fasta_writer[genome_name] = new fastqWriter(file_name);
        }
    }

    /*
        3- Processing reads
    */

    gzFile fp_1, fp_2;
    kseq_t* kseq_1, * kseq_2;
    fp_1 = gzopen(R1_file.c_str(), "r");
    fp_2 = gzopen(R2_file.c_str(), "r");

    kseq_1 = kseq_init(fp_1);
    kseq_2 = kseq_init(fp_2);

    int Reads_chunks_counter = 0;

    cout << "Partitioning the reads" << endl;
    if (write_fastq) {
        cout << "Writing the reads to fasta files" << endl;
    }
    else {
        cout << "Only stats" << endl;
    }

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

        double coverage_threshold = 0.1;
        double abundance_threshold = 2.0;
        uint32_t total_kmers = kmers_matches.size();
        
        std::pair<string, std::vector<uint32_t>> category;

        category = sensitive_classify_and_match_read_kmers(kmers_matches, stats);


        if (write_fastq) {
            if (category.first == "unique") {
                assert(category.second.size() == 1);
                fasta_writer[id_to_genome[category.second[0]]]->write(kseq_1, kseq_2);
            }
        }
    }

    // Write stats
    stats.print_json_to_file(base_filename + ".json");

}
