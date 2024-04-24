#ifndef RNANUE_DETECT_HPP
#define RNANUE_DETECT_HPP

// Standard
#include <iostream>
#include <bitset>

// Boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

// Classes
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
};

#endif //RNANUE_DETECT_HPP
