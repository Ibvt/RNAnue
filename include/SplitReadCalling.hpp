#ifndef RNANUE_DETECT_HPP
#define RNANUE_DETECT_HPP

//boost
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>

// openMP
#include <omp.h>


#include <iostream>
#include <fstream>
#include <regex>
#include <list>
#include <bitset>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <ctime>

// seqan3
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/alignment_file/all.hpp>


#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>

// filters
#include "ScoringMatrix.hpp"
#include "Traceback.hpp"

namespace po = boost::program_options;
namespace pt = boost::property_tree;
namespace fs = boost::filesystem;

using seqan3::operator""_tag;
using seqan3::operator""_cigar_op;
using seqan3::operator""_dna5;
using seqan3::operator""_dna15;
using seqan3::operator""_dna4;
using seqan3::get;

// overload struct to
template <> struct seqan3::sam_tag_type<"XX"_tag> { using type = int32_t; };
template <> struct seqan3::sam_tag_type<"XY"_tag> { using type = int32_t; };
template <> struct seqan3::sam_tag_type<"XJ"_tag> { using type = int32_t; };
template <> struct seqan3::sam_tag_type<"XH"_tag> { using type = int32_t; };

// create custom tags: complementarity (FC), energy hybridization (FE)
//template <> struct seqan3::sam_tag_type<"FC"_tag> { using type = std::string; };
//template <> struct seqan3::sam_tag_type<"FE"_tag> { using type = std::string; };

//template <> struct seqan3::sam_tag_type<"XC"_tag> { using type = std::string; };
//template <> struct seqan3::sam_tag_type<"XE"_tag> { using type = std::string; };
template <> struct seqan3::sam_tag_type<"XM"_tag> { using type = int32_t; }; // matches in alignment
template <> struct seqan3::sam_tag_type<"XL"_tag> { using type = int32_t; }; // length of alignment
template <> struct seqan3::sam_tag_type<"XN"_tag> { using type = int32_t; };

template <> struct seqan3::sam_tag_type<"XS"_tag> { using type = std::string; };


typedef std::pair<uint32_t,uint32_t> ReadPos;
typedef std::pair<uint64_t,uint64_t> GenomePos;
typedef std::vector<seqan3::cigar> CigarSplt;


// introduce record_type
using types = seqan3::type_list<
        std::string,
        seqan3::sam_flag,
        std::optional<int32_t>,
        std::optional<int32_t>,
        std::optional<uint8_t>,
        std::vector<seqan3::cigar>,
        std::span<seqan3::dna5>,
        seqan3::sam_tag_dictionary>;

using types_as_ids = seqan3::fields<
        seqan3::field::id,
        seqan3::field::flag,
        seqan3::field::ref_id,
        seqan3::field::ref_offset,
        seqan3::field::mapq,
        seqan3::field::cigar,
        seqan3::field::seq,
        seqan3::field::tags>;

using SamRecord = seqan3::record<types, types_as_ids>;
using ComplResult = std::tuple<int, int, double, double, std::vector<char>, std::vector<char>>;


typedef std::vector<
        std::tuple<
                std::string,
                seqan3::sam_flag,
                std::optional<int32_t>,
                std::optional<int32_t>,
                std::vector<seqan3::cigar>,
                seqan3::dna5_vector,
                seqan3::sam_tag_dictionary>> Splts;


class SplitReadCalling {
private:
    po::variables_map params;
    std::vector<std::tuple<std::string>> stats;

    int readscount;
    int alignedcount;
    int splitscount;
    int msplitscount;
    int nsurvivedcount;

    std::map<std::string, int> frequency;
    std::map<std::string, std::vector<std::pair<std::pair<int,int>,std::string>>> features;


public:
    SplitReadCalling();
    SplitReadCalling(po::variables_map params);
    ~SplitReadCalling();

    // iterate through reads
    void iterate(std::string matched, std::string splits, std::string multsplits);
    void process(auto &splitrecords, auto &splitsfile, auto &multsplitsfile);
    void distribute(auto &subrecords, auto &splits, auto &msplits);

    void filterSegments(auto &splitrecord, std::optional<int32_t> &refOffset,
                        std::vector<seqan3::cigar> &cigar, std::span<seqan3::dna5> &seq,
                        seqan3::sam_tag_dictionary &tags, std::vector<SamRecord> &curated);

    void addFilterToSamRecord(SamRecord &rec, std::pair<float,float> filters);
    void addComplementarityToSamRecord(SamRecord &rec1, SamRecord &rec2, TracebackResult &res);
    void addHybEnergyToSamRecord(SamRecord &rec1, SamRecord &rec2, double &hyb);

    void writeSamFile(auto &samfile, std::vector<std::pair<SamRecord,SamRecord>> &splits );

    TracebackResult complementarity(std::span<seqan3::dna5> &seq1, std::span<seqan3::dna5> &seq2);
    double hybridize(std::span<seqan3::dna5> &seq1, std::span<seqan3::dna5> &seq2);
    void createDir(fs::path path);

    void progress(std::ostream& out);

    int countSamEntries(std::string file, std::string command);
    std::vector<std::vector<fs::path>> splitInputFile(std::string matched, std::string splits, int entries);

    std::string addSuffix(std::string _file, std::string _suffix, std::vector<std::string> _keys);
    void start(pt::ptree sample);
};



/*
#include <iostream>
// Boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

// Standard
#include "IBPTree.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

class SplitReadCalling {
    public:
        SplitReadCalling(po::variables_map params);
        ~SplitReadCalling();

        void start(pt::ptree sample);
        void iterate(std::string& matched, std::string& splits, std::string& multsplits);

    private:
        po::variables_map params;
        IBPTree features;

};*/

#endif //RNANUE_DETECT_HPP
