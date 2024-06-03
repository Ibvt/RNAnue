#ifndef RNANUE_ANALYSIS_HPP
#define RNANUE_ANALYSIS_HPP

// Standard
#include <iostream>

// Boost
#include <boost/program_options.hpp>

namespace po = boost::program_options;

class Analysis {
    public:
        Analysis(po::variables_map params);
        ~Analysis();

        void start();
};

#endif //RNANUE_ANALYSIS_HPP