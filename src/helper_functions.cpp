#include <cstdint>
#include <chrono>
#include <glob.h>
#include <sstream>
#include <string>
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
#include <helper_functions.hpp>
#include <map>


using namespace std;


vector<file_info> file_fetcher(const std::string& pattern) {

    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if (return_value != 0) {
        globfree(&glob_result);
        std::stringstream ss;
        ss << "glob() failed with return_value " << return_value << std::endl;
        throw std::runtime_error(ss.str());
    }

    std::vector<file_info> filenames;
    for (size_t i = 0; i < glob_result.gl_pathc; ++i) {

        file_info finfo;

        finfo.path = std::string(glob_result.gl_pathv[i]);
        size_t _lastindex = finfo.path.find_last_of(".");
        finfo.prefix = finfo.path.substr(0, _lastindex);
        finfo.base_name = finfo.path.substr(finfo.path.find_last_of("/\\") + 1);

        filenames.push_back(finfo);
    }
    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}


void load_all_bins(string bins_dir, BINS_PHMAP* bin_to_hashes, int cores) {
    auto all_files = file_fetcher(bins_dir + "/*.bin");
    cout << "Fetched " << all_files.size() << " genome(s) .." << endl;
    bin_to_hashes->reserve(all_files.size());
#pragma omp parallel num_threads(cores)
    {
#pragma omp for
        for (int x = 0; x < all_files.size(); x++) {
            const auto& file = all_files[x];
            size_t lastindex = file.path.find_last_of(".");

            std::string::size_type idx = file.path.rfind('.');
            std::string extension = "";

            if (idx != std::string::npos) extension = file.path.substr(idx + 1);
            if (extension != "bin") {
                cerr << "skipping " << file.path << " does not have extension .bin" << endl;
                continue;
            }

            phmap::flat_hash_map<uint64_t, uint32_t> bin_hashes;
            phmap::BinaryInputArchive ar_in(file.path.c_str());
            bin_hashes.phmap_load(ar_in);

            auto bin_size = bin_hashes.size();
            auto _basename = file.base_name.substr(0, file.base_name.find_last_of("."));

            bin_to_hashes->try_emplace_l(_basename,
                [](BINS_PHMAP::value_type& v) {},
                bin_hashes
            );


#pragma omp critical
            cout << "Loaded " << _basename << endl;
        }
    }
}

string create_dir(string output_file, int serial) {
    int dir_err;
    string new_name = "";

    if (!serial) {
        dir_err = mkdir(output_file.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        new_name = output_file;
    }
    else {
        new_name = output_file + "_v." + to_string(serial);
        dir_err = mkdir(new_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    if (-1 == dir_err) return create_dir(output_file, ++serial);

    return new_name;
}


inline string time_diff(std::chrono::high_resolution_clock::time_point& t1) {
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



/*
The classify_and_match_read_kmers function takes the following parameters:
- genome_ids: A vector containing genome IDs matched with read kmers.
- coverage_threshold: The minimum percentage of read kmers that should match a genome to consider it significant.
  This helps to filter out genomes with a low number of matching kmers, which may be due to random chance or sequencing errors.
- ratio_threshold: The minimum ratio of read kmers of the dominant genome (the genome with the most matching kmers)
  compared to the second most common genome. This helps to ensure that the dominant genome is significantly more
  abundant than the others, increasing confidence in the classification.
*/
std::pair<string, std::vector<uint32_t>> classify_and_match_read_kmers(const std::vector<uint32_t>& genome_ids, Stats& stats, double coverage_threshold, double ratio_threshold) {
    stats.increment_reads();
    if (genome_ids.empty()) {
        stats.increment_unmapped();
        return { "unmapped", {} };
    }

    flat_hash_map<uint32_t, uint32_t> genome_id_count;
    for (const auto& id : genome_ids) {
        genome_id_count[id]++;
    }

    auto max_it = std::max_element(genome_id_count.begin(), genome_id_count.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });

    double coverage = static_cast<double>(max_it->second) / genome_ids.size();
    // cout << "coverage: " << coverage << " threshold: " << coverage_threshold << "\n";
    if (coverage < coverage_threshold) {
        stats.increment_unmapped();
        return { "unmapped", {} };
    }

    std::vector<uint32_t> matching_genome_ids{max_it->first};
    genome_id_count.erase(max_it);

    for (const auto& [genome_id, count] : genome_id_count) {
        double ratio = static_cast<double>(max_it->second) / count;
        // cout << "ratio(" << genome_id << "): " << ratio << " threshold: " << ratio_threshold << "\n";
        if (ratio < ratio_threshold) {
            matching_genome_ids.push_back(genome_id);
        }
    }

    if (matching_genome_ids.size() > 1) {
        for (const auto& id : matching_genome_ids) {
            stats.increment_ambiguous(id);
        }
        return { "ambiguous", matching_genome_ids };
    }
    else {
        stats.increment_unique(matching_genome_ids[0]);
        return { "unique", matching_genome_ids };
    }
}


// Coverage on total kmers
std::pair<string, std::vector<uint32_t>> classify_and_match_read_kmers(const std::vector<uint32_t>& genome_ids, uint32_t total_kmers, Stats& stats, double coverage_threshold, double ratio_threshold) {
    stats.increment_reads();
    if (genome_ids.empty()) {
        stats.increment_unmapped();
        return { "unmapped", {} };
    }

    flat_hash_map<uint32_t, uint32_t> genome_id_count;
    for (const auto& id : genome_ids) {
        genome_id_count[id]++;
    }

    auto max_it = std::max_element(genome_id_count.begin(), genome_id_count.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });

    double coverage = static_cast<double>(max_it->second) / total_kmers;
    // cout << "coverage: " << coverage << " threshold: " << coverage_threshold << "\n";
    if (coverage < coverage_threshold) {
        stats.increment_unmapped();
        return { "unmapped", {} };
    }

    std::vector<uint32_t> matching_genome_ids{max_it->first};
    genome_id_count.erase(max_it);

    for (const auto& [genome_id, count] : genome_id_count) {
        double ratio = static_cast<double>(max_it->second) / count;
        // cout << "ratio(" << genome_id << "): " << ratio << " threshold: " << ratio_threshold << "\n";
        if (ratio < ratio_threshold) {
            matching_genome_ids.push_back(genome_id);
        }
    }

    if (matching_genome_ids.size() > 1) {
        for (const auto& id : matching_genome_ids) {
            stats.increment_ambiguous(id);
        }
        return { "ambiguous", matching_genome_ids };
    }
    else {
        stats.increment_unique(matching_genome_ids[0]);
        return { "unique", matching_genome_ids };
    }
}


// Without stats
std::pair<string, std::vector<uint32_t>> classify_and_match_read_kmers(const std::vector<uint32_t>& genome_ids, double coverage_threshold, double ratio_threshold) {
    if (genome_ids.empty()) {
        return { "unmapped", {} };
    }

    flat_hash_map<uint32_t, uint32_t> genome_id_count;
    for (const auto& id : genome_ids) {
        genome_id_count[id]++;
    }

    auto max_it = std::max_element(genome_id_count.begin(), genome_id_count.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });

    double coverage = static_cast<double>(max_it->second) / genome_ids.size();
    // cout << "coverage: " << coverage << " threshold: " << coverage_threshold << "\n";
    if (coverage < coverage_threshold) {
        return { "unmapped", {} };
    }

    std::vector<uint32_t> matching_genome_ids{max_it->first};
    genome_id_count.erase(max_it);

    for (const auto& [genome_id, count] : genome_id_count) {
        double ratio = static_cast<double>(max_it->second) / count;
        // cout << "ratio(" << genome_id << "): " << ratio << " threshold: " << ratio_threshold << "\n";
        if (ratio < ratio_threshold) {
            matching_genome_ids.push_back(genome_id);
        }
    }

    if (matching_genome_ids.size() > 1) {
        return { "ambiguous", matching_genome_ids };
    }
    else {
        return { "unique", matching_genome_ids };
    }
}

// Super Sensitive
std::pair<string, std::vector<uint32_t>> sensitive_classify_and_match_read_kmers(const std::vector<uint32_t>& genome_ids, Stats& stats) {
    stats.increment_reads();
    if (genome_ids.empty()) {
        stats.increment_unmapped();
        return { "unmapped", {} };
    }

    flat_hash_map<uint32_t, uint32_t> genome_id_count;
    for (const auto& id : genome_ids) {
        genome_id_count[id]++;
    }

    auto max_it = std::max_element(genome_id_count.begin(), genome_id_count.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });

    std::vector<uint32_t> matching_genome_ids{max_it->first};
    genome_id_count.erase(max_it);

    for (const auto& [genome_id, count] : genome_id_count) {
        matching_genome_ids.push_back(genome_id);
    }

    if (matching_genome_ids.size() > 1) {
        for (const auto& id : matching_genome_ids) {
            stats.increment_ambiguous(id);
        }
        return { "ambiguous", matching_genome_ids };
    }
    else {
        stats.increment_unique(matching_genome_ids[0]);
        return { "unique", matching_genome_ids };
    }
}



// STATS


void Stats::increment_unmapped() {
    this->unmatched++;
}

void Stats::increment_ambiguous(uint32_t genome_id) {
    this->stats[this->get_genome_name(genome_id)]["ambiguous"]++;
}

void Stats::increment_unique(uint32_t genome_id) {
    this->stats[this->get_genome_name(genome_id)]["unique"]++;
}

void Stats::set_id_to_genome(flat_hash_map<int, string>& external_id_to_genome) {
    for (const auto& [id, genome_name] : external_id_to_genome) {
        this->id_to_genome[id] = genome_name;
        this->stats[genome_name] = {
            {"ambiguous", 0},
            {"unique", 0}
        };
    }
}

string Stats::get_genome_name(uint32_t id) {
    return id_to_genome[id];
}

void Stats::increment_reads(){
    this->reads++;
}


void Stats::print_json_to_file(string output_file) {
    ofstream out(output_file);
    out << "{" << endl;
    out << "\t\"total_reads\": " << reads << "," << endl;
    out << "\t\"unmapped\": " << unmatched << "," << endl;
    out << "\t\"unmapped%\": " << 100*((double)unmatched/reads) << "," << endl;
    out << "\t\"mapped\": {" << endl;
    for (const auto& [genome_name, genome_stats] : stats) {
        out << "\t\t\"" << genome_name << "\": {" << endl;
        for (const auto& [matching_class, count] : genome_stats) {
            out << "\t\t\t\"" << matching_class << "\": " << count << "," << endl;
            out << "\t\t\t\"" << matching_class << "%\": " << 100*((double)count/reads) << "," << endl;
        }
        out << "\t\t}," << endl;
    }

    out << "\t}" << endl;
    out << "}" << endl;
    out.close();
}
