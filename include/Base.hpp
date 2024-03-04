//
// Created by Richard Albin Schaefer on 1/23/24.
//

#ifndef RNANUE_BASE_HPP
#define RNANUE_BASE_HPP

#include <boost/program_options.hpp>
#include "Utility.hpp"
#include "Data.hpp"

namespace po = boost::program_options;

class Base {
    public:
        Base(po::variables_map params);
        ~Base();
        Data data;

    private:
        po::variables_map params;
//        Data data;
};

#endif //RNANUE_BASE_HPP
