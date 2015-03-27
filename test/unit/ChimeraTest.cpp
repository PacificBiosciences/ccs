#define BOOST_TEST_MODULE EndToEnd_Test

#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>
#include <zlib.h>

//#include "../../include/laamm/fastaparser/kseq.h"
//#include "../../include/laamm/fastaparser/FastaEntry.h"
//#include "../../include/laamm/chimera/ChimeraLabeler.h"
//#include "../../include/laamm/chimera/ChimeraLabel.h"
//#include "../../include/laamm/chimera/FastaSizePair.h"

#include <boost/test/unit_test.hpp>
#include <seqan/seq_io.h>

//using namespace PacBio::Chimera;

const std::string FILENAME = "../test/unit/data/test.fasta";

typedef std::pair<seqan::CharString, seqan::Dna5String> FastaEntry;

/**
 * @brief Parses a fasta file.
 * 
 * @param  File name.
 * @return Vector of FastaEntry
 */
std::vector<FastaEntry> ParseSingleFastaFile(const std::string&);

/**
 * @brief Create a list of adapter containing the SMRT-bell adapter.
 * @return shared_ptr with adapter
 */
std::shared_ptr<std::vector<FastaEntry>> CreateAdapterlist();

/**
 * @brief Create a vector of strings containing char-delimited substrings
 * @return vector of strings
*/
std::vector<std::string> split(const std::string &s, char delim);

BOOST_AUTO_TEST_CASE(fasta)
{
    // Get consensus sequences
    auto seqs = ParseSingleFastaFile(FILENAME);

    // Create instance of the ChimeraLabeler class, parameterized
    //     with the default minimum chimera score. 
    //ChimeraLabeler chimeraLabeler(1.0);
    //auto pairList = std::make_shared<std::vector<FastaSizePair>>();

    // For each consensus read
    //for (const auto& seq : seqs)
    //{
        //std::vector<std::string> idParts = split(seq.id, '_');
        //std::string numReadsString = idParts[3].substr(8);
        //int numReads = std::stoi(numReadsString);

        //std::cerr << seq.id << " " << numReads << std::endl;

        //FastaSizePair pair = std::make_pair(seq, numReads);
        //pairList->push_back(std::move(pair));
    //}
    //std::cerr << pairList->size() << std::endl;

    // Sort the pairs by size, in descending order
    //std::sort(pairList->begin(), pairList->end(), sortFastaSizePairs());
    //std::reverse(pairList->begin(), pairList->end());

    //chimeraLabeler.Label(pairList);
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<FastaEntry> ParseSingleFastaFile(std::string filePath)
{
    std::vector<FastaEntry> output;

    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;
    seqan::StringSet<seqan::CharString> quals;

    seqan::SeqFileIn seqFileIn(filePath);

    // Reads all remaining records.
    readRecords(ids, seqs, quals, seqFileIn);

    for (size_t i = 0; i < length(seqs); ++i)
    {
        FastaEntry fasta = std::make_pair(ids[i], seqs[i]);
        output.push_back(fasta);
    }

    return output;
}
