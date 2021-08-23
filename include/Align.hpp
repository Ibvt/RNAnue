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
#include <utility>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <seqan3/io/alignment_file/input.hpp>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>

#include <bitset>

namespace pt = boost::property_tree;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using seqan3::operator""_tag;
using seqan3::operator""_cigar_op;
using seqan3::operator""_dna5;
using seqan3::operator""_dna15;
using seqan3::operator""_dna4;
using seqan3::get;


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
        //seqan3::sam_tag_dictionary>> Splts;


//        seqan3::dna5_vector>> Splts;
       // seqan3::sam_tag_dictionary>
        //> Splts;*/

class Align {
    private:
        po::variables_map params;
        std::string index;


    public:
        // constructor
        Align(po::variables_map params);
        Align();

        void buildIndex();
        void alignReads(std::string query, std::string matched);
        void detSplits(std::string matched, std::string splits);

        //double complementarity(std::string rna1, std::string rna2);

        double complementarity(seqan3::dna5_vector rna1, seqan3::dna5_vector rna2);
        double hybridize(seqan3::dna5_vector rna1, seqan3::dna5_vector rna2);
		
        void processSplits(auto &splitrecords, auto &splitsfile);
        //std::vector<seqan3::dna5> spanToVec(std::span<seqan3::dna5,-1> seq);

        seqan3::dna5 string2dna5(std::string rna); 

        void constructIndex();
        void start(pt::ptree sample);
};
