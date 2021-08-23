#include "Base.hpp"

Base::Base(po::variables_map _params, std::string _subcall):
    params(_params),
    data(_params) {

        if(_params["subcall"].as<std::string>() == "preproc") {
            std::cout << "preproc called within Base" << std::endl;
            
            //
            if(_params["preproc"].as<std::bitset<1>>() != std::bitset<1>("1")) {
                std::cout << "'preproc' has not been set correctly - check parameters" << std::endl;
                exit(EXIT_FAILURE);
            }
            data.preproc();
        }

        if(_params["subcall"].as<std::string>() == "align") {
            std::cout << "align called within Base" << '\n';
            data.align();
        }

        if(_params["subcall"].as<std::string>() == "detect") {
            data.splitReadCalling();
        }

        if(_params["subcall"].as<std::string>() == "clustering") {
            std::cout << "cluster the split reads" << '\n';
			data.clustering();
        }

        if(_params["subcall"].as<std::string>() == "analysis") {
            std::cout << "analysis" << std::endl;
            data.analysis();
        }

        if(_params["subcall"].as<std::string>() == "complete") {
            std::cout << "complete analysis" << std::endl;
            data.preproc();
            data.align();
            data.splitReadCalling();
            data.analysis();
        }
}

//
template <typename Callable>
void Base::testBind(Callable f) {
    std::cout << "asda" << std::endl;

    pt::ptree bla;
    f();
}

