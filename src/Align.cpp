//
// Created by Richard Albin Schaefer on 2/12/24.
//

#include "Align.hpp"

Align::Align(po::variables_map params) : params(params) {
    // build index (if necessary)
    buildIndex();

}
Align::~Align() {}

void Align::buildIndex() {
    // get path of specified reference genome
    std::string ref = params["dbref"].as<std::string>();

    // determine path to store the index (within output directory)
    fs::path outDir = fs::path(params["outdir"].as<std::string>());
    fs::path idx = outDir / fs::path(ref).replace_extension(".idx").filename();

    if(fs::exists(idx)) {
        std::cout << helper::getTime() << "The reference index " << idx << " already exists\n";
    } else {
        // create index
        std::cout << helper::getTime() << "Create index for " << ref << "\n";
        std::string idxCall = "segemehl.x -x " + idx.string() + " -d " + ref;
        const char* idxCallChar = idxCall.c_str();
        system(idxCallChar);
    }
    index = idx.string();
}

//
void Align::alignReads(std::string query, std::string mate, std::string matched) {
    std::cout << helper::getTime() << "Start Alignment\n";
    std::string align = "segemehl.x";
    align += " -S "; // split mode
    align += " -b "; // output in bam format
    align += " -A " + std::to_string(params["accuracy"].as<int>());
    align += " -U " + std::to_string(params["minfragsco"].as<int>());
    align += " -W " + std::to_string(params["minsplicecov"].as<int>());
    align += " -Z " + std::to_string(params["minfraglen"].as<int>());
    align += " -t " + std::to_string(params["threads"].as<int>());
    align += " -m " + std::to_string(params["minlen"].as<int>());
    align += " -i " + index;
    align += " -d " + params["dbref"].as<std::string>();
    align += " -q " + query;
    if(mate.empty()) {
        align += " -p " + mate;
    }
    align += " -o " + matched;

    const char* call = align.c_str();
    system(call);
}


void Align::start(pt::ptree sample) {
    pt::ptree input = sample.get_child("input");
    pt::ptree output = sample.get_child("output");

    // align the reads (regular)
    std::string inFwd = input.get<std::string>("forward");
    std::string outMatched = output.get<std::string>("matched");

    alignReads(inFwd, "", outMatched);

    if(params["readtype"].as<std::string>() == "PE") {
        if(params["unprd"].as<std::bitset<1>>() == std::bitset<1>("1")) {
            std::string inR1only = input.get<std::string>("R1only");
            std::string outR1only = output.get<std::string>("matched_R1only");
            alignReads(inR1only, "", outR1only);

            std::string inR2only = input.get<std::string>("R2only");
            std::string outR2only = output.get<std::string>("matched_R2only");
            alignReads(inR2only, "", outR2only);
        }

        if(params["unmrg"].as<std::bitset<1>>() == std::bitset<1>("1")) {
            std::string inR1unmrg = input.get<std::string>("R1unmerged");
            std::string inR2unmrg = input.get<std::string>("R2unmerged");
            std::string outunmrg = output.get<std::string>("matched_unmerged");
            alignReads(inR1unmrg, inR2unmrg, outunmrg);
        }
    }
}


