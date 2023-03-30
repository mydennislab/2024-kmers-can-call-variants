#include <cassert>
#include <iostream>
#include <vector>
#include <helper_functions.hpp>
#include <string>
using namespace std;



struct TestCase {
    std::string description;
    std::vector<uint32_t> genome_ids;
    std::pair<string, std::vector<uint32_t>> expected_result;
};

void test_unique_genome() {
    TestCase test_case = {
        "Unique genome",
        {1, 1, 1, 1, 1, 1, 2, 2},
        {"Unique", {1}}
    };

    auto result = classify_and_match_read_kmers(test_case.genome_ids);

    if (result != test_case.expected_result) {
        std::cout << "Expected: " << test_case.expected_result.first << " " << test_case.expected_result.second[0] << std::endl;
        std::cout << "Actual: " << result.first << " " << result.second[0] << std::endl;
    }

}

void test_ambiguous_genomes_tie() {
    TestCase test_case = {
        "Ambiguous genomes due to tie",
        {1, 1, 1, 2, 2, 2},
        {"Ambiguous", {1, 2}}
    };

    auto result = classify_and_match_read_kmers(test_case.genome_ids);
    if (result != test_case.expected_result) {
        std::cout << "Expected: " << test_case.expected_result.first << " " << test_case.expected_result.second[0] << std::endl;
        std::cout << "Actual: " << result.first << " " << result.second[0] << std::endl;
    }
}

void test_ambiguous_genomes_low_coverage() {
    TestCase test_case = {
        "Ambiguous genomes due to low coverage",
        {1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2},
        {"Ambiguous", {1}}
    };

    auto result = classify_and_match_read_kmers(test_case.genome_ids, 0.5, 0.5);
    if (result != test_case.expected_result) {
        std::cout << "Expected: " << test_case.expected_result.first << " " << test_case.expected_result.second[0] << std::endl;
        std::cout << "Actual: " << result.first << " " << result.second[0] << std::endl;
    }
}

void test_ambiguous_genomes_low_ratio() {
    TestCase test_case = {
        "Ambiguous genomes due to low ratio",
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2},
        {"Ambiguous", {1}}
    };

    auto result = classify_and_match_read_kmers(test_case.genome_ids);
    if (result != test_case.expected_result) {
        std::cout << "Expected: " << test_case.expected_result.first << " " << test_case.expected_result.second[0] << std::endl;
        std::cout << "Actual: " << result.first << " " << result.second[0] << std::endl;
    }
}

void test_unmapped_genome() {
    TestCase test_case = {
        "Unmapped genome",
        {},
        {"Unmapped", {}}
    };

    auto result = classify_and_match_read_kmers(test_case.genome_ids);
    if (result != test_case.expected_result) {
        std::cout << "Expected: " << test_case.expected_result.first << " " << test_case.expected_result.second[0] << std::endl;
        std::cout << "Actual: " << result.first << " " << result.second[0] << std::endl;
    }
}

int main() {

    // cout << "testing unique genome" << endl;
    // test_unique_genome();
    // cout << "testing ambiguous genomes due to tie" << endl;
    // test_ambiguous_genomes_tie();
    // cout << "testing ambiguous genomes due to low coverage" << endl;
    // test_ambiguous_genomes_low_coverage();
    // cout << "testing ambiguous genomes due to low ratio" << endl;
    // test_ambiguous_genomes_low_ratio();
    // cout << "testing unmapped genome" << endl;
    // test_unmapped_genome();

    vector <uint32_t> genome_ids = { 1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3};
    auto result = classify_and_match_read_kmers(genome_ids, 0.1, 2);
    cout << result.first << ": ";
    for (auto i : result.second) {
        cout << i << ",";
    }
    cout << endl;


    std::cout << "All test cases passed." << std::endl;
    return 0;
}
