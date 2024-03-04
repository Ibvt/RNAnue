#include "SeqRickshaw.hpp"

SeqRickshaw::SeqRickshaw(po::variables_map params) {
    this->params = params;

    // get basic options for trimming
    int minlen = params["minlen"].as<int>();
    modetrm = params["modetrm"].as<int>();
    int qual = params["quality"].as<int>();
    threads = params["threads"].as<int>();

    // output
    std::string outDir = params["outdir"].as<std::string>();
    fs::path outDirSub = fs::path(outDir) / fs::path("preproc");

    // get adapter files
    std::string adpt5File = params["adpt5"].as<std::string>();
    std::string adpt3File = params["adpt3"].as<std::string>();

    try {
        switch(modetrm) {
            case 0: adpt5Tables = StateTransition("5\'-adapter", adpt5File);
                    break;
            case 1: {
                adpt3Tables = StateTransition("3\'-adapter", adpt3File);
                fs::path outDirSub3Adpt = outDirSub / fs::path("table_3adpt.txt");
                std::ofstream lookup3AdptOut(outDirSub3Adpt.string(), std::ios::trunc);
                adpt3Tables.writeStateTransitionTable(lookup3AdptOut);
                break;
            }

            case 2:
                adpt5Tables = StateTransition("5\'-adapter", adpt5File);
                adpt3Tables = StateTransition("3\'-adapter", adpt3File);
        }
    }
    catch (seqan3::file_open_error& err) {
        std::cout << err.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }
}

SeqRickshaw::~SeqRickshaw() {
}


void SeqRickshaw::start(pt::ptree sample) {

    // get the data objects from the sample
    pt::ptree input = sample.get_child("input");
    pt::ptree output = sample.get_child("output");

    fs::path inFwd = fs::path(input.get<std::string>("forward"));
    fs::path outFwd = fs::path(output.get<std::string>("forward"));

    Preproc pp = Preproc(params, adpt5Tables, adpt3Tables); // create preprocessor
    if(params["readtype"].as<std::string>() == "SE") {
        pp.processing(inFwd, outFwd); // process the reads (single-end)
    } else {
        // reverse read input/output
        fs::path inRev = fs::path(input.get<std::string>("reverse"));

        // filtering reads that are only in R1 or R2
        fs::path outR1only = fs::path(output.get<std::string>("R1only"));
        fs::path outR2only = fs::path(output.get<std::string>("R2only"));

        // unmerged reads
        fs::path outR1unmerged = fs::path(output.get<std::string>("R1unmerged"));
        fs::path outR2unmerged = fs::path(output.get<std::string>("R2unmerged"));

        // process the reads (paired-end)
        pp.processing(inFwd, outFwd, inRev, outR1only, outR2only, outR1unmerged, outR2unmerged);
    }
}
