#include <iostream>
#include <string>
#include <bitset>
#include <fstream>

// Boost
#include <boost/program_options.hpp>

// Seqan
#include <seqan3/core/debug_stream.hpp>

#include "Base.hpp"
#include "Closing.hpp"
#include "Config.hpp"
#include "Utility.hpp"

namespace po = boost::program_options;

void showVersion(std::ostream& _str) {
    _str << "RNAnue v" << RNAnue_VERSION_MAJOR;
    _str << "." << RNAnue_VERSION_MINOR << ".";
    _str << RNAnue_VERSION_PATCH << " - ";
    _str << "Detect RNA-RNA interactions ";
    _str << "from Direct-Duplex-Detection (DDD) data.";
    _str << std::endl;
}

void makePathAbs(po::variables_map& params, std::string param, fs::path& configFileDir) {
    if(params.count(param)) {
        fs::path filePath = fs::path(params[param].as<std::string>());
        if(filePath.empty()) { return; } // if file path is empty, do nothing
        if(!filePath.is_absolute()) {
            fs::path newPath = configFileDir / filePath;
            params.erase(param); // remove param from variables_map
            // add corrected path to variables_map
            params.insert(std::make_pair(param, po::variable_value((configFileDir / filePath).string(), false)));
        }
    }
}

// correct the paths to absolute paths (if they are relative)
void correctPaths(po::variables_map& params, fs::path& configFileDir) {
    std::vector<std::string> paramsToCheck = {"ctrls", "trtms", "outdir", "adpt3", "adpt5", "dbref", "features"};
    for(auto& param : paramsToCheck) {
        makePathAbs(params, param, configFileDir);
    }
}

int main(int argc, char* argv[]) {
    try {
        std::string readType;
        std::string configFile;

        po::options_description general("General");
        general.add_options()
                ("readtype,r", po::value<std::string>(&readType)->default_value("SE"),
                 "single-end (=SE) or paired-end (=PE) reads")
                ("trtms,t", po::value<std::string>()->default_value(""),
                 "folder containing the raw reads of the treatments including replicates located within subfolders (condition)")
                ("ctrls,s", po::value<std::string>()->default_value(""),
                 "folder containing the the raw reads of the controls including replicates located within subfolders (condition)")
                ("outdir,o", po::value<std::string>()->default_value(""),
                 "(output) folder in which the results are stored")
                ("threads,p", po::value<int>()->default_value(1),
                 "the number of threads")
                ("quality,q", po::value<int>()->default_value(20),
                 "lower limit for the quality (Phred Quality Score) of the reads")
                ("mapquality", po::value<int>()->default_value(20),
                 "lower limit for the quality (Phred Quality Score) of the alignments")
                ("minlen,l", po::value<int>()->default_value(15),
                 "minimum length of the reads")
                ("splicing", po::value<std::bitset<1>>()->default_value(0),
                 "splicing events are considered in the detection of split reads")
                ("features, f", po::value<std::string>()->default_value(""),
                        "annotation/features in .GFF3 format")
                ("overwrite", "overwrite existing results")
                ;

        po::options_description preproc("Preprocessing");
        preproc.add_options()
                ("preproc", po::value<std::bitset<1>>()->default_value(0),
                 "include preprocessing of the raw reads in the workflow of RNAnue")
                ("modetrm", po::value<int>()->default_value(1),
                 "mode of the trimming: only 5' (=0) and 3' (=1) or both (=2)")
                ("adpt5", po::value<std::string>()->default_value(""),
                 "file of the adapter sequences to be removed from the 5' end (.fasta)")
                ("adpt3", po::value<std::string>()->default_value(""),
                 "file of the adapter sequences to be removed from the 3' end (.fasta)")
                ("wtrim", po::value<std::bitset<1>>()->default_value(0),
                 "on whether (=1) or not (=0) to include window trimming")
                ("mmrate", po::value<double>()->default_value(0.1),
                 "rate of mismatched allowed when aligning adapter pattern to sequence")
                ("wsize", po::value<int>()->default_value(3),
                 "windows size to trim from 3' end")
                ("minovlps", po::value<int>()->default_value(5),
                 "minimal overlap to merge paired-end reads")
                ;

        po::options_description alignment("Alignment");
        alignment.add_options()
                ("dbref", po::value<std::string>(),
                 "reference genome (.fasta)")
                ("accuracy", po::value<int>()->default_value(90),
                 "minimum percentage of read matches")
                ("minfragsco", po::value<int>()->default_value(18),
                 "minimum score of a spliced fragment")
                ("minfraglen", po::value<int>()->default_value(20),
                 "minimum length of a spliced fragment")
                ("minsplicecov", po::value<int>()->default_value(80),
                        "minimum coverage for spliced transcripts")
                ("mapq", po::value<int>()->default_value(10),
                        "minimum mapping quality")
                ("exclclipping", po::value<std::bitset<1>>()->default_value(0),
                        "exclude soft clipping from the alignments")
                ("unprd", po::value<std::bitset<1>>()->default_value(0),
                 "only for paired-end reads: include unpaired reads")
                ("unmrg", po::value<std::bitset<1>>()->default_value(0),
                 "only for paired-end reads: include unmerged reads")
                ;

        po::options_description detect("Split Read Calling");
        detect.add_options()
                ("cmplmin", po::value<double>()->default_value(0.0), "complementarity cutoff for split reads")
                ("sitelenratio", po::value<double>()->default_value(0.1),
                 "aligned portion of the read (complementarity)")
                ("nrgmax", po::value<double>()->default_value(0), "hybridization energy cutoff for split reads")
                ;

        po::options_description clustering("Clustering");
        clustering.add_options()
                ("clust", po::value<std::bitset<1>>()->default_value(1),
                 "include clustering of the split reads in the workflow of RNAnue")
                ("clustdist", po::value<int>()->default_value(0), "minimum distance between clusters")
                ;

        po::options_description output("Output");
        output.add_options()
                ("stats", po::value<std::bitset<1>>()->default_value(0),
                 "whether (=1) or not (=0) to (additionally) create statistics of the libraries")
                ("outcnt", po::value<std::bitset<1>>()->default_value(0),
                 "whether (=1) or not (=0) to (additionally) save results as count table for DEA")
                ("outjgf", po::value<std::bitset<1>>()->default_value(0),
                 "whether (=1) or not (=0) to (additionally) save results as JSON Graph FILE (JGF)")
                ;

        po::options_description other("Other");
        other.add_options()
                ("version,v", "display the version number")
                ("help,h", "display this help message")
                ("config,c", po::value<std::string>(&configFile),
                 "configuration file that contains the parameters")
                ;

        po::options_description subcall("Subcall");
        subcall.add_options()
                ("subcall", po::value<std::string>(), "preproc, detect, alignment, clustering, analysis, complete")
                ;

        po::options_description cmdlineOptions;
        cmdlineOptions
                .add(general)
                .add(preproc)
                .add(alignment)
                .add(detect)
                .add(clustering)
                .add(output)
                .add(other)
                .add(subcall);

        po::options_description configFileOptions;
        configFileOptions
                .add(general)
                .add(preproc)
                .add(alignment)
                .add(detect)
                .add(clustering)
                .add(output);

        // translate all positional options into subcall options
        po::positional_options_description p;
        p.add("subcall", -1);

        po::variables_map params;
        store(po::command_line_parser(argc, argv).options(cmdlineOptions).positional(p).run(), params);
        notify(params);

        Closing cl;

        if(params.count("help")) {
            std::cout << cmdlineOptions << "\n";
            cl.printQuote(std::cout);
        }

        // check if subcall is empty
        if(params["subcall"].empty()) {
            std::cout << helper::getTime() << "Please provide a subcall\n";
            exit(EXIT_FAILURE);
        }
        showVersion(std::cout);

        // include parameters from the configfile
        std::ifstream ifs{configFile};

        if(params.count("version")) {
            showVersion(std::cout);
            cl.printQuote(std::cout);
            return 0;
        }

        if(!configFile.empty()) { // check that a config file has been provided
            std::ifstream ifs{configFile};
            if(!ifs) {
                std::cout << helper::getTime() << " Configuration file " << configFile << " could not be opened!\n";
                return 0;
            } else{
                // update the parameters from configfile
                po::store(po::parse_config_file(ifs, configFileOptions), params);
                notify(params);
            }
        }

        // correct the paths (if they are relative)
        fs::path configFileDir = fs::path(configFile).parent_path();
        correctPaths(params, configFileDir);

        Base base(params);

    } catch(po::error& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
    }
}
