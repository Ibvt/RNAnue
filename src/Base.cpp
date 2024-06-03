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
                    //std::cout << "cluster the split reads" << '\n';
                    //data.clustering();
                } else{
                    if(subcall == "analysis") {
                        //std::cout << "analysis" << std::endl;
                        //data.analysis();
                    } else {
                        if(subcall == "complete") {
                            data.preproc();
                            data.detect();
                            data.align();



                            //std::cout << "complete analysis" << std::endl;
                            /*
                            data.preproc();
                            data.align();
                            data.splitReadCalling();
                            data.analysis();*/
                        } else {
                            //std::cout << "subcall: " << _params["subcall"].as<std::string>() << " invalid!" << std::endl;
                            //exit(EXIT_FAILURE);
                        }
                    }
                }
            }
        }
    }
}

Base::~Base() {}


