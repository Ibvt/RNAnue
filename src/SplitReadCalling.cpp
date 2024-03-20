#include "SplitReadCalling.hpp"

//
SplitReadCalling::SplitReadCalling(po::variables_map params) {
    this->params = params;
    // create IBPTree (if splicing need to be considered)
    if(params["splicing"].as<std::bitset<1>>() ==  std::bitset<1>("1")) {
        // create IBPT
        this->features = IBPTree(params, 3);
    }
}

SplitReadCalling::~SplitReadCalling() {}


void SplitReadCalling::start(pt::ptree sample) {
    // input
    pt::ptree input = sample.get_child("input");
    std::string matched = input.get<std::string>("matched");

    // get other usefule variables
    int threads = params["threads"].as<int>();

    // output
    pt::ptree output = sample.get_child("output");
    std::string splits = output.get<std::string>("splits");
    std::string multsplits = output.get<std::string>("multsplits");



}






