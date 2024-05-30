#include "SplitReadCalling.hpp"

//
SplitReadCalling::SplitReadCalling(po::variables_map params) {
    this->params = params;
    // create IBPTree (if splicing need to be considered)
    if(params["splicing"].as<std::bitset<1>>() ==  std::bitset<1>("1")) {
        // order of the IBPTree
        int k = 7; // can be later extracted from config file
        // create IBPT
        this->features = IBPTree(params, k);
    }

    // initialize stats
    this->stats = Stats();



    const char* seq = "GGGAAAUCC";

}

SplitReadCalling::~SplitReadCalling() {}

void SplitReadCalling::iterate(std::string &matched, std::string &splits, std::string &multsplits) {
    // input BAM record
    seqan3::sam_file_input fin{matched};

    std::string currentQNAME = "";
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records;
    for(auto && record: fin) {

        // ignore reads if unmapped and do not suffice length requirements
        if(static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped)) { continue; }
        if(static_cast<int>(record.sequence().size()) < params["minlen"].as<int>()) { continue; }

        std::string QNAME = record.id();
        if((currentQNAME != "") && (currentQNAME != QNAME)) {

        }



    }


}


void SplitReadCalling::start(pt::ptree sample, pt::ptree condition) {
    // input
    pt::ptree input = sample.get_child("input");
    std::string matched = input.get<std::string>("matched");

    // print pt:ptree to see the structure
    std::string str = condition.data();

    std::cout << str << std::endl;



    // get other usefule variables
    int threads = params["threads"].as<int>();

    // output
    pt::ptree output = sample.get_child("output");
    std::string single = output.get<std::string>("single");
    std::string splits = output.get<std::string>("splits");
    std::string multsplits = output.get<std::string>("multsplits");

    /*
    std::cout << single << std::endl;
    std::cout << splits << std::endl;
    std::cout << multsplits << std::endl;
     */


    iterate(matched, splits, multsplits);
}






