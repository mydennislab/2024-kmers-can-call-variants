#include "string"
#include "kDataFrame.hpp"
#include <stdexcept>
#include "tuple"
#include "algorithms.hpp"
#include <sys/stat.h>
#include "colored_kDataFrame.hpp"
#include <vector>
#include <map>
#include <fstream>

using namespace std;

void eraseSubStr(std::string & mainStr, const std::string & toErase)
{
    // Search for the substring in string
    size_t pos = mainStr.find(toErase);
    if (pos != std::string::npos)
    {
        // If found then erase it from string
        mainStr.erase(pos, toErase.length());
    }
}

uint64_t kf_size(kDataFrame * KF){

    uint64_t size = 0;
    auto it = KF->begin();
    while(it != KF->end()){
        size++;
        it++;
    }
    return size;
}

std::string base_name(std::string const & path)
{
  return path.substr(path.find_last_of("/\\") + 1);
}

int main(int argc, char ** argv){

    // Loading genome KF
    string genome_prefix = argv[1];
    eraseSubStr(genome_prefix, ".mqf");
    string genome_basename = base_name(genome_prefix);
    cerr << "Loading genome KF ... ";
    auto genomeKF = kDataFrame::load(genome_prefix);
    cerr << "[DONE]" << endl;

    // Loading sample KF
    string samples_file = argv[2];
    eraseSubStr(samples_file, ".mqf");
    cout << "Loading " << samples_file << endl;
    string sample_basename = base_name(samples_file);
    cerr << "Loading Sample KF ... ";
    auto sampleKF = kDataFrame::load(samples_file);
    cerr << "[DONE]" << endl;

    cerr << "Infering sample size ..";
    uint64_t sample_kmers = kf_size(sampleKF);
    cerr << " [" << sample_kmers << "]" << endl;

    // Getting Intersection
    cerr << "calculating common kmers ("<< sample_basename << " & " << genome_prefix <<") ... ";
    auto commonKmersKF = kProcessor::kFrameIntersect({sampleKF, genomeKF});
    cerr << "[DONE]" << endl;
    uint64_t commonKmers = kf_size(commonKmersKF);
    double percentage = 100 * ((double)commonKmers / (double)sample_kmers);


    // Writing results
    cout << genome_basename << '\t' << sample_basename << '\t' << commonKmers << '\t' << sample_kmers << 't' << percentage << endl;
    cerr << "Writing results to " << "contamStats_" + genome_basename + "-" + sample_basename + ".tsv" << " ..." << endl;

    ofstream fs;
    string output_file_name = "contamStats_" + genome_basename + "-" + sample_basename + ".tsv";
    fs.open(output_file_name);
    fs << genome_basename << '\t' << sample_basename << '\t' << commonKmers << '\t' << sample_kmers << '\t' << percentage << endl;
    fs.close();
}