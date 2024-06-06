#ifndef RNANUE_ANALYSIS_HPP
#define RNANUE_ANALYSIS_HPP

// Standard
#include <iostream>

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
    private:
        po::variables_map params;
        IBPTree features;
};

#endif //RNANUE_ANALYSIS_HPP