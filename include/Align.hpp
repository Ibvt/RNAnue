// boost
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <regex>

#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/cigar/cigar_op.hpp>

#include <bitset>

namespace pt = boost::property_tree;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using seqan3::operator""_tag;
using seqan3::operator""_cigar_op;
using seqan3::operator""_dna5;
using seqan3::get;

// overload struct to 
template <> struct seqan3::sam_tag_type<"XX"_tag> { using type = int32_t; };
template <> struct seqan3::sam_tag_type<"XY"_tag> { using type = int32_t; };
template <> struct seqan3::sam_tag_type<"XJ"_tag> { using type = int32_t; };
template <> struct seqan3::sam_tag_type<"XH"_tag> { using type = int32_t; };

typedef std::pair<uint32_t,uint32_t> ReadPos;
typedef std::pair<uint64_t,uint64_t> GenomePos;
typedef std::vector<seqan3::cigar> CigarSplt;


typedef std::vector<
	std::tuple<
		std::string,
		seqan3::sam_flag, 
		std::optional<int32_t>,
		std::optional<int32_t>,
        std::vector<seqan3::cigar>,
        seqan3::dna5_vector,
        seqan3::sam_tag_dictionary>> Splts;

class Align {
    private:
        po::variables_map params;
        std::string index;


    public:
        // constructor
        Align(po::variables_map params);
        Align();

        void buildIndex();
        void alignReads(std::string query, std::string matched, std::string splits);
        void detSplits(std::string matched, std::string splits);

        double hybridize(std::string rna1, std::string rna2);
        double complementarity(std::string rna1, std::string rna2);
		
        void processSplits(auto &splitrecords, auto &splitsfile);
        std::vector<seqan3::dna5> spanToVec(std::span<seqan3::dna5,-1> seq);

        void constructIndex();
        void start(pt::ptree sample);
};

