#ifndef DATA_HPP
#define DATA_HPP

#include <iostream>
#include <fstream>
#include <limits>
#include <deque>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>


#include "SeqRickshaw.hpp"
#include "Align.hpp"
#include "SplitReadCalling.hpp"
#include "Clustering.hpp"
#include "Analysis.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

using GroupsPath = std::map<std::string, fs::path>; // path to ctrls, trtms

template <class DataType>
using GroupsMap = std::map<std::string, std::pair<std::string, DataType>>;

typedef std::vector<fs::path> PathVector;

class Data {
    /*
     * the type of the data:
     * 
    */
    private:
        po::variables_map params;
        pt::ptree dataStructure;

    public:
        Data();
        Data(po::variables_map _params);
        
        // getter & setter
        //
        pt::ptree getDataStructure();

        //
        void preprocDataPrep();
        void alignDataPrep();
        void detectDataPrep();
        void clusteringDataPrep();
        void analysisDataPrep();

        //  
        GroupsPath retrieveGroupsPath(fs::path _ctrls, fs::path _trtms);
        void retrieveData(GroupsPath _groupsPaths);
        pt::ptree retrieveGroup(std::string _group, fs::path _condition);
        pt::ptree retrieveOutput(fs::path _outConditionDir, pt::ptree _input);

        // iterates through the data structure and excutes the subcall
        
        template <typename Callable>
        void callInAndOut(Callable f);

        // helper methods
        // determine content of directory and sort it (return as vector)
        PathVector sortDirContent(fs::path _path);
		// filter content of directory to only include files containing search string
		PathVector filterDirContent(PathVector vec, std::string sestr);
        std::string addSuffix(std::string _file, std::string _suffix, std::vector<std::string> _keys);


        // test stuff
        template <typename Callable>
        void bla(Callable f);

        void preproc();
        void align();
        void splitReadCalling();
        void clustering();
        void analysis();



};

#endif // DATA_HPP
