// Standard
#include <iostream>
#include <bitset>

// Class
#include "Base.hpp"

Base::Base(po::variables_map params) : params(params), paramsVal(params), data(params) {
    std::string subcall = params["subcall"].as<std::string>();

    // preprocessing
    if(subcall == "preproc") {
        // preproc has not been actively switch on - 'preproc' parameter - but preproc has been called (subcall)
        if(params["preproc"].as<std::bitset<1>>() != std::bitset<1>("1")) {
            std::cout << helper::getTime() << " 'preproc' has not been set correctly - check parameters" << std::endl;
            exit(EXIT_FAILURE);
        }
        data.preproc();
    } else {
        if(subcall == "align") {
            data.align();
        } else {
            if(subcall == "detect") {
                data.detect();
            } else {
                if(subcall == "clustering") {
                    data.clustering();
                } else{
                    if(subcall == "analysis") {
                        data.analysis();
                    } else {
                        if(subcall == "complete") {
                            data.preproc();
                            data.detect();
                            data.align();
                            if(params["clust"].as<std::bitset<1>>() == std::bitset<1>("1")) {
                                data.clustering();
                            }
                            data.analysis();
                        } else {
                            std::cout << "subcall: " << params["subcall"].as<std::string>() << " invalid!" << std::endl;
                            exit(EXIT_FAILURE);
                        }
                    }
                }
            }
        }
    }
}

Base::~Base() {}


