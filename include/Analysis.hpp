#ifndef RNANUE_ANALYSIS_HPP
#define RNANUE_ANALYSIS_HPP

// Standard
#include <iostream>
#include <map>
#include <vector>

// Boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

// Class
#include "IBPTree.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

class Analysis {
    public:
        Analysis(po::variables_map params);
        ~Analysis();

        void start(pt::ptree sample);
        void iterate(std::string& single, std::string& splits);

    private:
        po::variables_map params;
        IBPTree features;
        std::map<std::pair<std::string, std::string>, double> pdf;
};

#endif //RNANUE_ANALYSIS_HPP