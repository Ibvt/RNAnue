#ifndef RNANUE_DATA_HPP
#define RNANUE_DATA_HPP

// Standard
#include <typeinfo>
#include <deque>

// Boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

// SeqAn3
#include "Utility.hpp"
#include "SeqRickshaw.hpp"
#include "Align.hpp"
#include "SplitReadCalling.hpp"
#include "Clustering.hpp"
#include "Analysis.hpp"

// namespaces
namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;
namespace jp = boost::property_tree::json_parser;

// map that contains the paths to the groups (e.g., ctrls, trtms)
using GroupsPath = std::map<std::string, fs::path>;
using PathVector = std::vector<fs::path>;

class Data{
    public:
        Data(po::variables_map params);
        ~Data();

        // get Data
        GroupsPath getGroupsPath(fs::path& ctrls, fs::path& trtms);
        void getCondition(GroupsPath& groups, std::string subcall);
        pt::ptree getData(const std::string group, fs::path& condition, std::string& subcall);

        // get output data
        pt::ptree getOutputData(pt::ptree& input, fs::path& conditionOutDir, std::string& subcall);
        pt::ptree getPreprocOutputData(pt::ptree& input, fs::path& conditionOutDir);
        pt::ptree getAlignOutputData(pt::ptree& input, fs::path& conditionOutDir);
        pt::ptree getDetectOutputData(pt::ptree& input, fs::path& conditionOutDir);
        pt::ptree getAnalysisOutputData(pt::ptree& input, fs::path& conditionOutDir);

        //
        int getNumberElements(PathVector& vec, std::string& subcall);
        std::vector<std::string> getSampleKeys(std::string& subcall);

        // prep functions
        void preprocDataPrep();
        void alignDataPrep();
        void detectDataPrep();
        void clusteringDataPrep();
        void analysisDataPrep();

        //
        template <typename Callable>
        void callInAndOut(Callable f, std::string subcallStr);

        // callables
        void preproc();
        void align();
        void detect();
        void clustering();
        void analysis();

    private:
        po::variables_map params;
        pt::ptree dataStructure;

};


#endif //RNANUE_DATA_HPP
