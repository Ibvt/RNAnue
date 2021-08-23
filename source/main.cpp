#include <iostream>
#include <string>
#include <bitset>

#include "Base.hpp"
#include "Config.h"
#include "Closing.hpp"

//#include <boost/program_options.hpp>

namespace po = boost::program_options;

void showVersion(std::ostream& _str) {
    _str << "RNAnue v" << RNAnue_VERSION_MAJOR;
    _str << "." << RNAnue_VERSION_MINOR << ".";
    _str << RNAnue_VERSION_PATCH << " - ";
    _str << "Detect RNA-RNA interactions ";
    _str << "from Direct-Duplex-Detection (DDD) data."; 
    _str << std::endl;
}

// A helper function to simplify the main part.
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " ")); 
    return os;
}

int main(int argc, char* argv[]) {
    try{

        std::string readType;
        std::string configFile;

        po::options_description general("General");
        general.add_options()
		    ("readtype,r", po::value<std::string>(&readType)->default_value("SE"), 
                "single-end (=SE) or paired-end (=PE) reads")
            ("trtms,t", po::value<std::string>(), 
                "folder containing the raw reads of the treatments including replicates located within subfolders (condition)")
            ("ctrls,s", po::value<std::string>(), 
                "folder containing the the raw reads of the controls including replicates located within subfolders (condition)")
            ("outdir,o", po::value<std::string>(), 
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
            ("exclclipping", po::value<std::bitset<1>>()->default_value(0),
                "exclude soft clipping from the alignments")
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

        po::options_description analysis("Analysis");
        analysis.add_options()
            ("features", po::value<std::string>(), "annotation/features in .GFF3 format")
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
            ("config,c", po::value<std::string>(&configFile)->default_value("params.cfg"),
                "configuration file that contains the parameters")
        ;

        po::options_description subcall("Subcall");
        subcall.add_options()
            ("subcall", po::value<std::string>(), "preproc, detect, alignment, clustering, analysis")
        ;

        po::options_description cmdlineOptions;
        cmdlineOptions
            .add(general)
            .add(preproc)
            .add(alignment)
            .add(detect)
            .add(clustering)
            .add(analysis)
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
            .add(analysis)
            .add(output);

        // translate all positional options into subcall options 
        po::positional_options_description p;
        p.add("subcall", -1);

        po::variables_map vm;
        store(po::command_line_parser(argc, argv).options(cmdlineOptions).positional(p).run(),vm);
        notify(vm);


        // include parameters from the configfile if available
        std::ifstream ifs(configFile.c_str());
    
        if(!ifs) {
            std::cout << "configuration file" << configFile << "could not be opened!" << std::endl;
            return 0;
        } else {
            po::store(po::parse_config_file(ifs, configFileOptions),vm);
            notify(vm);
        }


        Closing cl; // class that handles the closing remarks
        
        if(vm.count("help")) {
            std::cout << cmdlineOptions << std::endl;
            cl.show(std::cout);
            return 0;
        }

        if(vm.count("version")) {
            showVersion(std::cout);
            cl.show(std::cout);
            return 0;
        }

        // start execution
        showVersion(std::cout);
        Base bs(vm, vm["subcall"].as<std::string>()); // controls all downstream processing

        /* 
        if(vm.count("subcall")) {
            if(vm["subcall"].as<std::string>() == "preproc") {
                std::cout << "*** start the preprocessing" << std::endl;
            }

            if(vm["subcall"].as<std::string>() == "alignment") {
                std::cout << "alignment" << std::endl;
            }
            
            if(vm["subcall"].as<std::string>() == "clustering") {
                std::cout << "clustering" << std::endl;
            }
            
            if(vm["subcall"].as<std::string>() == "detect") {
                std::cout << "detect" << std::endl;
            }
            
            if(vm["subcall"].as<std::string>() == "complete") {
                std::cout << "perform complete analysis" << std::endl;
            }
        }*/

        cl.show(std::cout);
       
    
    } catch(po::error& e) {
        std::cerr << e.what();
        return 0;
    }
}

