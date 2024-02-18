#include "string"
#include <stdexcept>
#include "tuple"
#include <sys/stat.h>
#include <zlib.h>
#include <cstdio>
#include "kseq/kseq.h"
#include <iostream>
#include <vector>
#include <cassert>

using namespace std;
using namespace phmap;
KSEQ_INIT(gzFile, gzread)


string create_dir(string output_file, int serial) {
    int dir_err;
    string new_name = "";

    if (!serial) {
        dir_err = mkdir(output_file.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        new_name = output_file;
    } else {
        new_name = output_file + "_v." + to_string(serial);
        dir_err = mkdir(new_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    if (-1 == dir_err) return create_dir(output_file, ++serial);

    return new_name;
}


class fastqWriter {

public:
    ofstream fileStream_r1, fileStream_r2;

    explicit fastqWriter(string &filename_prefix) {
        this->fileStream_r1.open(filename_prefix + "_R1.fastq");
        this->fileStream_r2.open(filename_prefix + "_R2.fastq");
    }

    void write(kseq_t *seq1, kseq_t *seq2) {

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

inline bool file_exists(const std::string &name) {
    struct stat buffer
            {
            };
    return (stat(name.c_str(), &buffer) == 0);
}

inline string time_diff(std::chrono::high_resolution_clock::time_point &t1) {
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    long hr = milli / 3600000;
    milli = milli - 3600000 * hr;
    long min = milli / 60000;
    milli = milli - 60000 * min;
    long sec = milli / 1000;
    milli = milli - 1000 * sec;
    string timeDiff;
    timeDiff.append(to_string(min));
    timeDiff.append(":");
    timeDiff.append(to_string(sec));
    timeDiff.append(":");
    timeDiff.append(to_string(milli));

    return timeDiff;
}

flat_hash_map<uint64_t, std::vector<uint32_t>> load_colors(string index_prefix) {
    flat_hash_map<uint64_t, std::vector<uint32_t> > colors;
    string inputFilename = index_prefix + "colors.intvectors";
    ifstream input(inputFilename);
    uint32_t size;
    input >> size;
    colors = flat_hash_map<uint64_t, std::vector<uint32_t> >(size);
    for (int i = 0; i < size; i++) {
        uint64_t color, colorSize;
        input >> color >> colorSize;
        uint32_t sampleID;
        colors[color] = std::vector<uint32_t>(colorSize);
        for (int j = 0; j < colorSize; j++) {
            input >> sampleID;
            colors[color][j] = sampleID;
        }
    }
    return colors;
}

tuple<string, vector<int>> score(vector<uint32_t> &genomes) {

    vector<int> sources;

    if (genomes.empty())
        return make_tuple("unmapped", sources);

    flat_hash_map<int, int> scores;
    flat_hash_map<int, int> reverse_scores;
    flat_hash_map<int, int> countFreq;
    vector<int> all_scores;

    for (const auto &genome: genomes) {
        scores[genome]++;
    }

    for (const auto &score: scores) {
        countFreq[score.second]++;
        all_scores.emplace_back(score.second);
        reverse_scores[score.second] = score.first;
    }

    auto max = std::max_element(all_scores.begin(), all_scores.end());

    if (countFreq[*max] == 1) {
        sources.emplace_back(reverse_scores[*max]);
        return make_tuple("unique", sources);
    }

    for (const auto &score: scores) {
        if (score.second == *max) {
            sources.emplace_back(score.first);
        }
    }

    return make_tuple("ambig", sources);;

}


int main(int argc, char **argv) {

    if (argc != 4) {
        cerr << "run ./fastq_partitioner <index_prefix> <in_1.fasta> <in_2.fasta> <op:chunk_size def:1000>" << endl;
        exit(1);
    }

    string index_prefix = argv[1];
    string R1_file = argv[2];
    string R2_file = argv[3];


    if (!file_exists(R1_file)) {
        throw std::runtime_error("Could not open the unitigs fasta file");
    }

    if (!file_exists(R2_file)) {
        throw std::runtime_error("Could not open the unitigs fasta file");
    }

    if (!file_exists(index_prefix + ".extra")) {
        throw std::runtime_error("Could not open kProcessor index file");
    }

    // kProcessor Index Loading
    std::cerr << "Loading kProcessor index ..." << std::endl;
    colored_kDataFrame *ckf = colored_kDataFrame::load(index_prefix);
    kDataFrame *kf = ckf->getkDataFrame();


    set<int> vec_singleColors;
    flat_hash_map<uint64_t, vector<uint32_t>> color_to_vecGroups;
    flat_hash_map<uint64_t, string> color_to_groupString;
    cerr << "Loading colors .." << endl;
    auto colorsIntVector = load_colors(index_prefix);
    for (auto const &color : colorsIntVector) {
        uint64_t color_id = color.first;
        auto all_group_ids = color.second;

        for (auto _grp_id : all_group_ids) {
            color_to_vecGroups[color_id].emplace_back(_grp_id);
            vec_singleColors.emplace(_grp_id);
        }
    }

    map<string, fastqWriter *> fasta_writer;

    cerr << "Creating fasta file handlers" << endl;
    for (auto &item : vec_singleColors) {
        string file_name = "genome_" + to_string(item) + "_readsPartition";
        fasta_writer[to_string(item)] = new fastqWriter(file_name);
    }

    string unmapped_file_name = "unmapped_partition.fa";
    fasta_writer["unmapped"] = new fastqWriter(unmapped_file_name);


    int kSize = 21;
    int hashing_mode = 1;


    gzFile fp_1, fp_2;
    kseq_t *kseq_1, *kseq_2;
    fp_1 = gzopen(R1_file.c_str(), "r");
    fp_2 = gzopen(R2_file.c_str(), "r");

    kseq_1 = kseq_init(fp_1);
    kseq_2 = kseq_init(fp_2);

    cout << "Processing started ..." << endl;

    for (int seqCounter = 0; kseq_read(kseq_1) >= 0 && kseq_read(kseq_2) >= 0; seqCounter++) {

        uint32_t seq_1_length = string(kseq_1->seq.s).size();
        uint32_t seq_2_length = string(kseq_2->seq.s).size();

        if (seq_1_length < kSize || seq_2_length < kSize) continue;

        std::string seq = string(kseq_1->seq.s) + string(kseq_2->seq.s);

        vector<uint32_t> kmers_matches;

        for (unsigned long i = 0; i < seq.size() - kSize + 1; i++) {
            uint64_t color = kf->getCount(seq.substr(i, kSize));
            for (const auto &genomeID : color_to_vecGroups[color]) {
                kmers_matches.emplace_back(genomeID);
            }
        }

        auto category = score(kmers_matches);
        if (get<0>(category) == "unmapped") {
            fasta_writer["unmapped"]->write(kseq_1, kseq_2);
        } else if (get<0>(category) == "unique") {
            assert(get<1>(category).size() == 1);
            fasta_writer[to_string(get<1>(category)[0])]->write(kseq_1, kseq_2);
        } else if (get<0>(category) == "ambig") {
            for (auto const &genomeID : get<1>(category)) {
                fasta_writer[to_string(genomeID)]->write(kseq_1, kseq_2);
            }
        }


    }

    for (auto f : fasta_writer)
        f.second->close();

    kseq_destroy(kseq_1);
    kseq_destroy(kseq_2);
    gzclose(fp_1);
    gzclose(fp_2);


    return 0;
}