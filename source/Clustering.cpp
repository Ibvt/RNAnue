#include <iostream>
#include "Clustering.hpp"

Clustering::Clustering(po::variables_map params) : 
    params(params) {
    	std::cout << "create clustering object" << std::endl;
}



void Clustering::start(pt::ptree sample) {
    pt::ptree input = sample.get_child("input");
    pt::ptree output = sample.get_child("input");

    std::string splits = input.get<std::string>("splits");


    


    //std::cout << "within start" << splits << std::endl;
}
