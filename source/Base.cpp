#include "Base.hpp"

Base::Base(po::variables_map _params, std::string _subcall):
    params(_params),
    data(_params) {
        std::cout << "create Base" << std::endl;

        if(_params["subcall"].as<std::string>() == "preproc") {
            std::cout << "preproc called within Base" << std::endl;
            data.preproc();
        }

        if(_params["subcall"].as<std::string>() == "align") {
            std::cout << "align called within Base" << '\n';
            data.align();
        }

        if(_params["subcall"].as<std::string>() == "clustering") {
            std::cout << "cluster the split reads" << '\n';
        }

        if(_params["subcall"].as<std::string>() == "complete") {
            std::cout << "complete analysis" << std::endl;
        }
}

//
template <typename Callable>
void Base::testBind(Callable f) {
    std::cout << "asda" << std::endl;

    pt::ptree bla;
    f();
}

