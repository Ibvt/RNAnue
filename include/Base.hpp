#ifndef RNANUE_BASE_HPP
#define RNANUE_BASE_HPP

// Boost
#include <boost/program_options.hpp>

// Class
#include "Data.hpp"
#include "ParameterValidation.hpp"
#include "Utility.hpp"

namespace po = boost::program_options;

class Base {
    public:
        Base(po::variables_map params);
        ~Base();

    private:
        po::variables_map params;
        ParameterValidation paramsVal;
        Data data;
};

#endif //RNANUE_BASE_HPP
