
#ifndef RNANUE_ANALYSIS_HPP
#define RNANUE_ANALYSIS_HPP

#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <
#include <boost/filesystem.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <string>
#include<math.h>

#include <iostream>
#include <fstream>
#include <bitset>
#include <vector>
#include <tuple>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/range/views/chunk.hpp>

#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>

namespace pt = boost::property_tree;
namespace po = boost::program_options;
namespace fs = boost::filesystem;


using seqan3::operator""_tag;

template <> struct seqan3::sam_tag_type<"XE"_tag> { using type = std::string; };
template <> struct seqan3::sam_tag_type<"XC"_tag> { using type = std::string; };

class Analysis {

private:
    po::variables_map params;
    std::map<std::string, std::vector<std::pair<std::pair<int,int>,std::string>>> features;
    std::map<std::string, int> frequency;
    std::map<std::pair<std::string, std::string>, double> pdf;

    std::vector<std::string> interPaths; // paths of the interactions file

public:
    Analysis();
    Analysis(po::variables_map params);

    std::string retrieveTagValue(std::string tags, std::string tagName, std::string oldValue);
    void createCountTable();
    float calc_pvalue(int x, int n, float p);

    void parseAnnotations();
    void start(pt::ptree sample);

};

#endif // ANALYSIS_HPP

#endif //RNANUE_ANALYSIS_HPP
