#ifndef RNANUE_PREPROC_HPP
#define RNANUE_PREPROC_HPP

// Standard
#include <iostream>
#include <mutex>
#include <future>
#include <bitset>
#include <numeric> // std::accumulate
#include <ranges>
#include <vector>

// Boost
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

// Seqan
#include <seqan3/io/views/async_input_buffer.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/record.hpp>
#include <seqan3/alphabet/views/all.hpp>
#include <seqan3/utility/views/convert.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

// Classes
#include "StateTransition.hpp"
#include "Utility.hpp"

namespace po = boost::program_options;
namespace pt = boost::property_tree;

class Preproc {
    public:
        Preproc(po::variables_map& params, StateTransition& adpt5Tables, StateTransition& adpt3Tables);
        ~Preproc();

        // handling for single-end and paired-end reads
        void processing(fs::path& in, fs::path& out); // handling for single-end reads
        void processing(fs::path& inFwd, fs::path& outFwd, fs::path& inRev, fs::path& outR1only,
                        fs::path& outR2only, fs::path& outR1unmerged, fs::path& outR2unmerged);

        // trimming
        std::pair<std::size_t,std::size_t> trimReads(dtp::DNAVector& seq);
        std::size_t boyermoore(dtp::DNAVector& read, dtp::StateTransitionTable& table,
                          dtp::DNAVector bases, const dtp::DNAVector& pat);
        std::size_t nibble(dtp::QualVector& qual, std::pair<std::size_t,std::size_t>& bnds); // window trimming

        // merging
        void merging(std::string& fwdSeqid, auto& fwdSeq, auto& fwdQual, std::string& revSeqid, auto& revSeq,
                     auto& revQual, auto& outFwd, auto& outR1only, auto& outR2only, auto& outR1unmrg, auto& outR2unmrg);
        std::pair<dtp::DNAVector, dtp::QualVector> mergeReads(dtp::DNASpan& fwdSeq, dtp::QualSpan& fwdQual,
                                                                   dtp::DNASpan& revSeq, dtp::QualSpan& revQual);
        std::string longestCommonSubseq(std::string& fwd, std::string& rev);

        // filtering
        bool filterReads(auto& qual);
        double calcAvgPhred(auto& qual); // calculate the average phred score

    private:
        po::variables_map params;
        StateTransition adpt5Tables;
        StateTransition adpt3Tables;
        int wsize;
        int quality;
        int minlen;
        int minovlps;
};


#endif //RNANUE_PREPROC_HPP
