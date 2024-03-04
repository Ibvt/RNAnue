//
// Created by Richard Albin Schaefer on 1/24/24.
//

#ifndef RNANUE_DATA_HPP
#define RNANUE_DATA_HPP

#include <typeinfo>
#include <deque>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include "Utility.hpp"
#include "SeqRickshaw.hpp"
#include "SplitReadCalling.hpp"
#include "Align.hpp"

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
        void getCondition(GroupsPath& groups);
        pt::ptree getData(const std::string group, fs::path& condition);

        // get output data
        pt::ptree getOutputData(pt::ptree& input, fs::path& conditionOutDir);
        pt::ptree getPreprocOutputData(pt::ptree& input, fs::path& conditionOutDir);
        pt::ptree getAlignOutputData(pt::ptree& input, fs::path& conditionOutDir);
        pt::ptree getDetectOutputData(pt::ptree& input, fs::path& conditionOutDir);

        //
        int getNumberElements(PathVector& vec);
        std::vector<std::string> getSampleKeys();

        // prep functions
        void preprocDataPrep();
        void alignDataPrep();
        void detectDataPrep();

        //
        template <typename Callable>
        void callInAndOut(Callable f);

        // callables
        void preproc();
        void align();
        void detect();

    private:
        po::variables_map params;
        pt::ptree dataStructure;

};


#endif //RNANUE_DATA_HPP
