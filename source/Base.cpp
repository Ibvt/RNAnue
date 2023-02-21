#include "Base.hpp"

Base::Base(po::variables_map _params, std::string _subcall):
    params(_params),
    data(_params) {

        // preprocessing
        if(_params["subcall"].as<std::string>() == "preproc") {
            std::cout << "preproc" << std::endl;
            // preproc has not been actively switch on - 'preproc' parameter
            if(_params["preproc"].as<std::bitset<1>>() != std::bitset<1>("1")) {
                std::cout << "'preproc' has not been set correctly - check parameters" << std::endl;
                exit(EXIT_FAILURE);
            }
            data.preproc();
        } else {
            if(_params["subcall"].as<std::string>() == "align") {
                std::cout << "align called within Base" << '\n';
                data.align();
            } else {
                if(_params["subcall"].as<std::string>() == "detect") {
                    data.splitReadCalling();
                } else {
                    if(_params["subcall"].as<std::string>() == "clustering") {
                        std::cout << "cluster the split reads" << '\n';
                        data.clustering();
                    } else{
                        if(_params["subcall"].as<std::string>() == "analysis") {
                            std::cout << "analysis" << std::endl;
                            data.analysis();
                        } else {
                            if(_params["subcall"].as<std::string>() == "complete") {
                                std::cout << "complete analysis" << std::endl;
                                data.preproc();
                                data.align();
                                data.splitReadCalling();
                                data.analysis();
                            } else {
                                std::cout << "subcall: " << _params["subcall"].as<std::string>() << " invalid!" << std::endl;
                                exit(EXIT_FAILURE);
                            }
                        }
                    }
                }
            }
        }
}

//
template <typename Callable>
void Base::testBind(Callable f) {
    std::cout << "asda" << std::endl;

    pt::ptree bla;
    f();
}

