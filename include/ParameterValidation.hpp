#ifndef RNANUE_PARAMETERVALIDATION_HPP
#define RNANUE_PARAMETERVALIDATION_HPP

// Standard
#include <iostream>

// Boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// Class
#include "Utility.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

class ParameterValidation {
    public:
        ParameterValidation(po::variables_map& params);
        ~ParameterValidation();

        // validation routines
        void validatePreproc();
        void validateDirs(std::string param);

    private:
        po::variables_map& params;
};

#endif //RNANUE_PARAMETERVALIDATION_HPP
