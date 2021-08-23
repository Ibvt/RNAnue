#include "Align.hpp"

Align::Align(po::variables_map params) : 
    params(params) {

    std::cout << "create align object" << std::endl;
    buildIndex();
}

void Align::alignReads(std::string query, std::string matched) {
    std::string align = "segemehl.x";
    align += " -S ";
    align += " -A " + std::to_string(params["accuracy"].as<int>()); 
    align += " -U " + std::to_string(params["minfragsco"].as<int>());
    align += " -W " + std::to_string(params["minsplicecov"].as<int>());
    align += " -Z " + std::to_string(params["minfraglen"].as<int>());
    align += " -t " + std::to_string(params["threads"].as<int>());
    align += " -m " + std::to_string(params["minlen"].as<int>());
    align += " -i " + index;
    align += " -d " + params["dbref"].as<std::string>();
    align += " -q " + query;
    align += " -o " + matched;

    const char* call = align.c_str();
    system(call);
    std::cout << align << std::endl;
}

void Align::buildIndex() {
    // retrieve path of reference genome
    std::string ref = params["dbref"].as<std::string>();

    fs::path outDir = fs::path(params["outdir"].as<std::string>());
    fs::path gen = outDir / fs::path(ref).replace_extension(".idx").filename();

    if(fs::exists(gen)) {
        std::cout << "segemehl index found on filesystem\n";
    } else {
        std::cout << "generate index " << "\n";
        std::string genIndex = "segemehl.x -x " + gen.string() + " -d " + ref;
        const char* call = genIndex.c_str();
        system(call);
    }
    index = gen.string();
}


//
seqan3::dna5 Align::string2dna5(std::string rna) {
    seqan3::dna5 seq{};
    for(auto &nt : rna) {
        seq.assign_char(nt);
    }
    return seq;
}


//
void Align::start(pt::ptree sample) {
    pt::ptree input = sample.get_child("input");
    pt::ptree output = sample.get_child("output");

    std::cout << "start align" << std::endl;

    if(params["readtype"].as<std::string>() == "SE") {
        std::string forward = input.get<std::string>("forward");
        std::string matched = output.get<std::string>("matched");
        alignReads(forward, matched);

    }
}
    
