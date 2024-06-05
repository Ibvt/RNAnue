#include "ParameterValidation.hpp"

ParameterValidation::ParameterValidation(po::variables_map& params) : params(params) {
    if(params["subcall"].as<std::string>().empty()) {
        std::cout << helper::getTime() << "Please provide a subcall\n";
        exit(EXIT_FAILURE);
    } else {
        // make sure a correct subcall has been provided
        std::string subcall = params["subcall"].as<std::string>();
        if(subcall == "preproc") {
            validatePreproc();
        } else if (subcall == "align") {

        } else if (subcall == "detect") {

        } else if (subcall == "clustering") {

        } else if (subcall == "analysis") {

        } else {
            std::cout << helper::getTime() << "Invalid subcall: " << subcall << "\n";
            exit(EXIT_FAILURE);
        }
    }

    if(params["subcall"].as<std::string>() == "preproc") {
        validatePreproc();
    }
}
ParameterValidation::~ParameterValidation() {}

void ParameterValidation::validatePreproc() {
    // validate input (trtms) and output (outdir) directories
    validateDirs("trtms");
    validateDirs("outdir");
}

void ParameterValidation::validateDirs(std::string param) {
    if (params[param].as<std::string>().empty()) {
        std::cout << helper::getTime() << "Please provide a value for --" << param << "\n";
        exit(EXIT_FAILURE);
    } else {
        if (param == "outdir") {
            // check if the output directory already exists
            if (fs::exists(params[param].as<std::string>())) {
                // if outdir already exists - check for overwrite flag
                /*
                if (!params.count("overwrite")) {
                    std::cout << helper::getTime() << "Specified output directory (--outdir ";
                    std::cout << params["outdir"].as<std::string>() << ") already exists\n";
                    std::cout << helper::getTime() << "Please provide a new directory or use the --overwrite flag\n";
                }*/
            } else {
                // check if the input directory exists
                if (!fs::exists(params[param].as<std::string>())) {
                    std::cout << helper::getTime() << "Specified output directory (-- " << param;
                    std::cout << params[param].as<std::string>() << ") does not exist\n";
                }
            }
        }
    }
}
