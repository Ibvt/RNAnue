#ifndef RNANUE_ALIGN_HPP
#define RNANUE_ALIGN_HPP

// Standard
#include <iostream>
#include <filesystem>
#include <bitset>

// Boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/process.hpp>

// Class
#include "Utility.hpp"

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

class Align {
    public:
        // constructor/destructor
        Align(po::variables_map params);
        ~Align();

        // alignment
        void start(pt::ptree sample);
        void buildIndex();
        void alignReads(std::string query, std::string mate, std::string matched);

    private:
        po::variables_map params;
        std::string index;
};

#endif //RNANUE_ALIGN_HPP
