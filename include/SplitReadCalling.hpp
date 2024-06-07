#ifndef RNANUE_DETECT_HPP
#define RNANUE_DETECT_HPP

// Standard
#include <iostream>
#include <bitset>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <thread>
#include <mutex>
#include <future>
#include <vector>
#include <iostream>
#include <sstream>
#include <regex>
#include <cstdlib>

// Boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

// HTSlib
#include <htslib/sam.h>
#include <htslib/hts.h>

// SeqAn3
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>

// Class
#include "IBPTree.hpp"
#include "Stats.hpp"
#include "FilterScores.hpp"
#include "ScoringMatrix.hpp"
#include "Traceback.hpp"

using types = seqan3::type_list<
        std::string,
        seqan3::sam_flag,
        std::optional<int32_t>,
        std::optional<int32_t>,
        dtp::CigarVector,
        dtp::DNASpan,
        dtp::QualSpan,
        seqan3::sam_tag_dictionary>;

using fields = seqan3::fields<
        seqan3::field::id,
        seqan3::field::flag,
        seqan3::field::ref_id,
        seqan3::field::ref_offset,
        seqan3::field::cigar,
        seqan3::field::seq,
        seqan3::field::qual,
        seqan3::field::tags>;
using SAMrecord = seqan3::sam_record<types, fields>;

// define tags
using seqan3::operator""_tag;
template <> struct seqan3::sam_tag_type<"XH"_tag> { using type = int32_t; }; // number of splits (SAM record)
template <> struct seqan3::sam_tag_type<"XX"_tag> { using type = int32_t; }; // begin of split
template <> struct seqan3::sam_tag_type<"XY"_tag> { using type = int32_t; }; // end of split
template <> struct seqan3::sam_tag_type<"XJ"_tag> { using type = int32_t; }; // number of splits (whole read)
template <> struct seqan3::sam_tag_type<"XN"_tag> { using type = int32_t; };
template <> struct seqan3::sam_tag_type<"XM"_tag> { using type = int32_t; }; // number of matches (complementarity)
template <> struct seqan3::sam_tag_type<"XL"_tag> { using type = int32_t; }; // length of alignment
template <> struct seqan3::sam_tag_type<"XC"_tag> { using type = float; }; // complementarity
template <> struct seqan3::sam_tag_type<"XR"_tag> { using type = float; }; // sitelenratio
template <> struct seqan3::sam_tag_type<"XA"_tag> { using type = std::string; }; // alignment
template <> struct seqan3::sam_tag_type<"XS"_tag> { using type = int32_t; }; // quality
template <> struct seqan3::sam_tag_type<"XE"_tag> { using type = float; }; // end of split


using namespace seqan3::literals;
using seqan3::get;

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

class SplitReadCalling {
    public:
        SplitReadCalling(po::variables_map params);
        ~SplitReadCalling();

        // main function
        void start(pt::ptree sample, pt::ptree condition);
        void sort(const std::string& inputFile, const std::string& outputFile); // sort BAM file


        void iterate(std::string& matched, std::string& single, std::string& splits, std::string& multsplits);
        template <typename T>
        void process(std::vector<T>& readrecords, auto& singleOut, auto& splitsOut, auto& multsplitsOut,
                     auto& singleOutMutex, auto& splitsOutMutex, auto& multsplitsOutMutex);
        void decide(std::map<int, std::vector<SAMrecord>>& putative, auto& splitsOut, auto& multsplitsOut,
                    auto& splitsOutMutex, auto& multsplitsOutMutex);
        bool filter(auto& sequence, uint32_t cigarmatch);
        bool matchSpliceSites(dtp::Interval& spliceSites, std::optional<uint32_t> refId);
        void storeSegments(auto& splitrecord, std::optional<int32_t> refPos, dtp::DNASpan& seq,
                           dtp::QualSpan& qual, dtp::CigarVector& cigar, seqan3::sam_tag_dictionary& tags,
                           std::vector<SAMrecord>& curated);

        // filters
        TracebackResult complementarity(dtp::DNASpan &seq1, dtp::DNASpan &seq2);
        double hybridization(dtp::DNASpan &seq1, dtp::DNASpan& seq2);

        // output
        void addComplementarityToSamRecord(SAMrecord &rec1, SAMrecord &rec2, TracebackResult &res);
        void addHybEnergyToSamRecord(SAMrecord &rec1, SAMrecord &rec2, double &hyb);
        void writeSAMrecordToBAM(auto& bamfile, std::vector<std::pair<SAMrecord, SAMrecord>>& records);


    private:
        po::variables_map params;
        IBPTree features;
        //Stats stats;
        std::shared_ptr<Stats> stats;
        std::string condition; // stores the current condition
        std::deque<std::string> refIds; // stores the reference ids
        FilterScores filterScores;
};

#endif //RNANUE_DETECT_HPP
