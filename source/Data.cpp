#include "Data.hpp"
#include "Helper.hpp"


//
void createDirectory(fs::path& _path, std::ostream& _str) {
    if(fs::exists(_path)) {
        _str << "\t" << _path << " already exists." << std::endl;
    } else {
        _str << "\t created " << _path << std::endl;
        fs::create_directory(_path);
    }
}

void createOutDir(fs::path& _path, std::string _subcall, std::ostream& _str) {
    _str << "*** create directories to store the results (specified via --outdir)" << std::endl;
    //_str << "\tcreate: " << _path.string() << std::endl;
    createDirectory(_path, _str);
    fs::path pathSub = _path / _subcall;
    createDirectory(pathSub, _str);
}


//
Data::Data() {}
Data::Data(po::variables_map _params) :
    params(_params) {
    // determine which subcall (e.g., preproc, align, clustering, analysis) was made 
    std::string subcall = _params["subcall"].as<std::string>();

    // create output directory
    fs::path resultsDir = fs::path(_params["outdir"].as<std::string>()); 
    createOutDir(resultsDir, subcall, std::cout);

    // preprocessing
    if(subcall == "preproc") { //
        preprocDataPrep(); // prepares the data files as property tree
    }

    if(subcall == "align") {
        alignDataPrep();
    }

    if(subcall == "detect") {
        detectDataPrep();
    }

	if(subcall == "clustering") {
		clusteringDataPrep();
	}

    if(subcall == "analysis") {
        analysisDataPrep();
    }


}

// getter & setter
pt::ptree Data::getDataStructure() {
    return this->dataStructure;
}

//
void Data::preprocDataPrep() {
    std::cout << "reads the data" << std::endl;
    // retrieve paths that contain all the reads
    fs::path ctrlsPath = fs::path(this->params["ctrls"].as<std::string>());
    fs::path trtmsPath = fs::path(this->params["trtms"].as<std::string>());

    /*
     * trtms -> "/Users/..../trtms"
     * (ctrls -> "/Users/.../crtls")
     * */
    GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath); // 
    retrieveData(group);
}


void Data::alignDataPrep() {
    std::cout << "reads the data to align" << '\n';

    if(params["preproc"].as<std::bitset<1>>() == std::bitset<1>(0)) {
        std::cout << "skipping preprocessing" << '\n';

        fs::path ctrlsPath = fs::path(params["ctrls"].as<std::string>());
        fs::path trtmsPath = fs::path(params["trtms"].as<std::string>());

        /*
         * trtms -> "/Users/..../trtms"
         * (ctrls -> "/Users/.../crtls")
         * */
        GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath); // 
        retrieveData(group);
    } else {
        if(params["preproc"].as<std::bitset<1>>() == std::bitset<1>(1)) {

            fs::path ctrlsPath = fs::path(params["outdir"].as<std::string>()) / "preproc/ctrls"; 
            fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / "preproc/trtms"; 
        
            GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath); // 
//            std::cout << "after groupspath" << std::endl;
            retrieveData(group);
        }
    }
}


void Data::detectDataPrep() {
    if(params["stats"].as<std::bitset<1>>() == 1) {
        fs::path statsfile = fs::path(params["outdir"].as<std::string>()) / "detectStat.txt";
        std::ofstream ofs;
        ofs.open(statsfile.string());
        ofs << "library\tmapped\tsplits\tmultisplits" << std::endl;
        ofs.close();
    }
    
    // determine path of the alignments results
    fs::path ctrlsPath = fs::path(params["outdir"].as<std::string>()) / "align/ctrls";
    if(params["ctrls"].as<std::string>() == "") {
        ctrlsPath = fs::path("");
    }
    fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / "align/trtms";
	
    GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath);
	retrieveData(group);
}

void Data::clusteringDataPrep() {
	std::cout << "reads the data to clustering" << "\n";
	
	fs::path ctrlsPath = fs::path(params["outdir"].as<std::string>()) / "detect/ctrls";
	if(params["ctrls"].as<std::string>() == "") {
		ctrlsPath = fs::path("");
	}
	fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / "detect/trtms";

	GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath);
	retrieveData(group);
}

void Data::analysisDataPrep() {
    std::cout << "reads the data for analysis" << "\n";

    fs::path ctrlsPath = fs::path(params["outdir"].as<std::string>()) / "detect/ctrls";
    if(params["ctrls"].as<std::string>() == "") {
        ctrlsPath = fs::path("");
    }
    fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / "detect/trtms";

    GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath);
    retrieveData(group);
}

//
GroupsPath Data::retrieveGroupsPath(fs::path _ctrls, fs::path _trtms) {
    GroupsPath groups;
    std::cout << helper::getTime() << " collect the datasets" << std::endl;

    if(_trtms != "") { // '--trtms' has been set in the config file or cmdline
        if(_ctrls != "") { // '--ctrls' has been set in the config or cmdline
            groups.insert(std::make_pair("ctrls", _ctrls));
        } else { // '--ctrls' has not been set in either config or cmdline - still continue
            std::cout << "### WARNING - this call runs RNAnue without control data" << std::endl;
            std::cout << "### WARNING - \"--ctrls\" has not been set" << std::endl;
        }
        groups.insert(std::make_pair("trtms", _trtms));
    } else { // abort - trtms have not been set in either config file or cmdline
        std::cerr << "### ERROR - RNAnue is aborted! No treatment data specified" << std::endl;
        exit(EXIT_FAILURE);
    }
    return groups;
}


// creates property of the files in the system
void Data::retrieveData(GroupsPath _groupsPath) {
    // ptree object containing the filesystem with data
    pt::ptree data, subcall, path, group, condition;

    // vector of fs::path that contains the conditions (as paths)
    PathVector conditionsVec;
    fs::path pathOut; // buffer output path

    //
    GroupsPath::iterator itGroups = _groupsPath.begin(); 
    for(itGroups;itGroups != _groupsPath.end();++itGroups) { // iterate through the groups (e.g., ctrls, trtms)
        // print paths to ctrls and trtms (e.g., trtms:     "/Users/.../trtms/)
        std::cout << "\t" << itGroups->first << ":\t" << itGroups->second << std::endl;
        if(fs::exists(itGroups->second)) { // check if the provided path exists...
            std::cout << "\t\t...has been found in the filesystem!" << std::endl;

            if(fs::is_directory(itGroups->second)) { // ... and check if its a directory (not a file)

                // retrieve content (conditions) within group (e.g., rpl-exp,..) as abs path
                PathVector conditionsVec = sortDirContent(itGroups->second); // sort

                PathVector::iterator itConditions = conditionsVec.begin();
                for(itConditions;itConditions != conditionsVec.end();++itConditions) {

                    if(fs::exists(*itConditions) && !fs::is_directory(*itConditions)) {
                        std::cout << "\t\tBut contains files (which will be ignored!)" << std::endl;
                        continue;
                    }
                    condition = retrieveGroup(itGroups->first,*itConditions);

                    // parent output directories
                    group.push_back(std::make_pair("", condition));
                }
                subcall.add_child(itGroups->first, group);   
                group.erase("");
            } else {
                std::cout << "### ERROR - " << itGroups->second << " is not a directory!" << std::endl;
                exit(EXIT_FAILURE);
            }
        } else {
            std::cout << "### ERROR - " << itGroups->second <<  "could not be found in the filesystem!" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    dataStructure.add_child(params["subcall"].as<std::string>(), subcall);


    std::cout << "*** write to output file" << std::endl;
    // write folder structure to output file
    fs::path ppPath = params["outdir"].as<std::string>() / fs::path("data.json");
    std::ofstream ppStr(ppPath.string());
    boost::property_tree::json_parser::write_json(ppStr, dataStructure);
}


// constructs property tree for a group (e.g., ctrls, trtms)
// parameters: e.g., trtms (_group), /Users/.../trtms/rpl_37C (_conditionPath)
pt::ptree Data::retrieveGroup(std::string _group, fs::path _conditionPath) {
    // define variables
    std::string subcall = params["subcall"].as<std::string>();

    int nrElements; // count of elements for each "dataset"  
    int smplCntr = 0;
    int elCntr;

    std::string conditionName = _conditionPath.stem().string();

    // define propery trees for the data
    pt::ptree condition, samples, sample, files, input, output;

	//std::cout << "_conditionPath " << _conditionPath << std::endl;
    
    std::vector<std::string> sampleKeys;
	
	// lists and sorts content in PathVector
    PathVector dataFiles = sortDirContent(_conditionPath);
	
    std::cout << "\t\t" << conditionName << std::endl;

    if(subcall == "preproc") {
        nrElements = (params["readtype"].as<std::string>() == "PE") ? 2 : 1;
        sampleKeys = {"forward", "reverse"}; 
    }
    
    if(subcall == "align") {
        if(params["preproc"].as<std::bitset<1>>() == std::bitset<1>(0)) {
            nrElements = 1;
            sampleKeys = {"forward"};
        } else {
            PathVector tmpDataFiles = filterDirContent(dataFiles, "preproc.fastq");

            nrElements = 1; // 
    //        sampleKeys = {"matched","splits"};
            sampleKeys = {"forward"};
            dataFiles = tmpDataFiles;
        }
    }

    if(subcall == "detect") {
        nrElements = 1;
        // remove files that do not contain keyword (matched)
        PathVector tmpDataFiles = filterDirContent(dataFiles, "matched.sam");
        sampleKeys = {"matched"};
        dataFiles = tmpDataFiles;
    }


	if(subcall == "clustering") {
		//std::cout << "being in retrieve group" << std::endl;
		nrElements = 1;

		// remove files that do not contain keyword (splits)
		PathVector tmpDataFiles = filterDirContent(dataFiles, "_splits.sam");
		sampleKeys = {"splits"};
		dataFiles = tmpDataFiles;
	}

    if(subcall == "analysis") {
        nrElements = 1;
		PathVector tmpDataFiles = filterDirContent(dataFiles, "_splits.sam");
		sampleKeys = {"splits"};
		dataFiles = tmpDataFiles;
    }

    /* path of the results up to the level of the conditions 
     * e.g., /results/preproc/trtms/
     * */
    fs::path outConditionDir = fs::path(params["outdir"].as<std::string>());
    outConditionDir /= fs::path(params["subcall"].as<std::string>());
    outConditionDir /= fs::path(_group);
    outConditionDir /= fs::path(_conditionPath.filename());

    //
    elCntr = 0;
    for(unsigned i=0;i<dataFiles.size();++i) {
        files.put(sampleKeys[i % nrElements],dataFiles[i].string());
        //
        if(elCntr == nrElements-1) { 
//            sample.put("group", _group);
 //           sample.put("condition", conditionName);
            // input and output
            sample.add_child("input", files);
            output = retrieveOutput(outConditionDir, sample);


            sample.add_child("output",output);
   
            samples.push_back(std::make_pair("",sample));
            sample.erase("input");
            sample.erase("output");
            for(unsigned j=0;j<sampleKeys.size();++j) {
                files.erase(sampleKeys[j]);
            }
            elCntr = 0;
        } else {
            ++elCntr;
        }
    }
    condition.put("condition", conditionName);
 //   condition.put("path", _condition.string());
    condition.add_child("samples",samples);
    
    return condition;
}


// switch path (up to level condition from input to Output
fs::path replacePath(fs::path _replacement, fs::path _original ) {
    fs::path nPath = _replacement / _original.filename();
    return nPath;
}

// create ptree of the output files
pt::ptree Data::retrieveOutput(fs::path _outConditionDir, pt::ptree _input) {
    pt::ptree output;
    // preprocessing
    if(params["subcall"].as<std::string>() == "preproc") {
        // replace input path with output path (results/...)
        fs::path fwd = fs::path(_input.get<std::string>("input.forward"));
        std::string forward = replacePath(_outConditionDir, fwd).string();
        
        // create output files using the input files with an additonal suffix
        std::string forwardOut = addSuffix(forward,"_preproc", {"_R1","fwd"});
        output.put("forward", forwardOut); // write output back to ptree

        // using paired-end reads results in additional files for unassembled reads
        if(params["readtype"].as<std::string>() == "PE") { 
            // replace 
            fs::path rev = fs::path(_input.get<std::string>("input.reverse"));
            std::string reverse = replacePath(_outConditionDir, rev).string();

            // create outfiles using the input files with an additional suffix
            std::string r1only = addSuffix(reverse,"_R1only", {"_R2","rev"});
            std::string r2only = addSuffix(reverse,"_R2only", {"_R2","rev"});

            std::string r1unmerged = addSuffix(reverse, "_R1unmerged", {"_R2","rev"});
            std::string r2unmerged = addSuffix(reverse, "_R2unmerged", {"_R2","rev"});

            output.put("R1only", r1only); // push both back to output ptree
            output.put("R2only", r2only);

            output.put("R1unmerged", r1unmerged);
            output.put("R2unmerged", r2unmerged);
        }

        
    }

    // alignment
    if(params["subcall"].as<std::string>() == "align") {
        std::cout << "retrieve output " << std::endl;

        //
        fs::path fwd = fs::path(_input.get<std::string>("input.forward"));
        
        fwd.replace_extension(".sam");
        std::string forward = replacePath(_outConditionDir, fwd).string();

        std::string matched = addSuffix(forward, "_matched", {});
 //       std::string splits = addSuffix(forward, "_splits", {});

        output.put("matched", matched);
//        output.put("splits", splits);

    }

    if(params["subcall"].as<std::string>() == "detect") {
//        output.put("splits","adads");

        fs::path matched = fs::path(_input.get<std::string>("input.matched"));
        std::string s = replacePath(_outConditionDir, matched).string();
        
        std::string splits = addSuffix(s, "splits", {"matched"});
        std::string multsplits = addSuffix(s, "multsplits", {"matched"});
        output.put("splits",splits);
        output.put("multsplits",multsplits);
    }

	// clustering
	if(params["subcall"].as<std::string>() == "clustering") {
        fs::path splits = fs::path(_input.get<std::string>("input.splits"));
		splits.replace_extension(".txt");
		std::string clu = replacePath(_outConditionDir, splits).string();
		std::string clusters = addSuffix(clu, "_clusters", {"_splits"});
        //output.put("clusters", clusters); # no individual outputs for each - instead one file
	}
    
    if(params["subcall"].as<std::string>() == "analysis") {
        fs::path splits = fs::path(_input.get<std::string>("input.splits"));
        splits.replace_extension(".txt");
        std::string its = replacePath(_outConditionDir, splits).string();
        std::string ints = addSuffix(its, "_interactions",{"_splits"});
        output.put("interactions", ints);
    }
    return output;
}

//
template <typename Callable>
void Data::callInAndOut(Callable f) {

    // retrieve paths for parameters
    fs::path outDir = fs::path(params["outdir"].as<std::string>());
    std::string subcallStr = params["subcall"].as<std::string>();
    fs::path outSubcallDir = outDir / fs::path(subcallStr);

    // create output directory (and subdirectory for subcall)
    createDirectory(outDir, std::cout); // create output
    createDirectory(outSubcallDir, std::cout); //  create subcall

    pt::ptree subcall = dataStructure.get_child(subcallStr);
    std::deque<std::string> groups = {"trtms"};
    if(subcall.size() > 1) {
        groups.push_front("ctrls");
    }
    
    pt::ptree conditions, samples, path;
    fs::path outGroupDir, outConditionDir; 

    for(unsigned i=0;i<groups.size();++i) {
        // create directory for groups (e.g., ctrls, trtms)
        outGroupDir = outSubcallDir / fs::path(groups[i]);

        // don't create subfolders for clustering && analysis
        //if(params["subcall"].as<std::string>() != "clustering") {
        createDirectory(outGroupDir, std::cout);
        //}

        conditions = subcall.get_child(groups[i]);
        BOOST_FOREACH(pt::ptree::value_type const &v, conditions.get_child("")) {
            pt::ptree condition = v.second;
            
            // create directory for condition (e.g., rpl_exp) 
            outConditionDir = outGroupDir / fs::path(condition.get<std::string>("condition"));
            createDirectory(outConditionDir, std::cout);
           
            samples = condition.get_child("samples");
            // iterate over samples
            BOOST_FOREACH(pt::ptree::value_type const &w, samples.get_child("")) { 
                pt::ptree sample = w.second;

                // call start function
                f(sample);
                //
                // f(path, group  sample)
            }
        }
    }
}

template <typename Callable>
void Data::bla(Callable f) {
    std::cout << "bla test" << std::endl;
    f();
}

void Data::preproc() {
    // create Object of Trimming
    SeqRickshaw srs(params);
    callInAndOut(std::bind(&SeqRickshaw::start, srs, std::placeholders::_1));
}

void Data::align() {
    // create Object of Alignment
    Align aln(params);
    callInAndOut(std::bind(&Align::start, aln, std::placeholders::_1));
}


void Data::splitReadCalling() {
    SplitReadCalling src(params);
    callInAndOut(std::bind(&SplitReadCalling::start, src, std::placeholders::_1));
}


void Data::clustering() {
    // create Object of CLustering
	Clustering clu(params);
	callInAndOut(std::bind(&Clustering::start, &clu, std::placeholders::_1));

    clu.sumup();
}

void Data::analysis() {
    Analysis anl(params);
    callInAndOut(std::bind(&Analysis::start, &anl, std::placeholders::_1));

    if(params["outcnt"].as<std::bitset<1>>() == std::bitset<1>(1)) {
        anl.createCountTable();
    }
}


/*
 * Helper Function
 */

// lists and sorts the content of a directory
PathVector Data::sortDirContent(fs::path _path) {
    PathVector content; // 
    copy(fs::directory_iterator(_path), fs::directory_iterator(), back_inserter(content));
    sort(content.begin(),content.end()); // sort the content 
    return content;
}


// filters content of directory
PathVector Data::filterDirContent(PathVector vec, std::string sestr) {
	PathVector content; // create new vector
	std::copy_if(vec.begin(), vec.end(), std::back_inserter(content), [&sestr] (fs::path x) {
		// only return if fs::path contains sestr
		return x.string().find(sestr) != std::string::npos; 
	});
	return content;
}



// adds suffix to filename (optional: 
std::string Data::addSuffix(std::string _file, std::string _suffix, std::vector<std::string> _keywords) {
    int keyPos, tmpKeyPos = -1; // buffer the positions of the keywords
    int dotPos = _file.find("."); // determine position of the dot
    keyPos = dotPos;
    if(!_keywords.empty()) {
        for(unsigned i=0;i<_keywords.size();++i) {
            tmpKeyPos = _file.find(_keywords[i]);
                if(tmpKeyPos != -1) { // key could be found
                    keyPos = tmpKeyPos;
                }
        }
    }
    std::string newFile = _file.substr(0,keyPos);
    newFile += _suffix;
    newFile += _file.substr(dotPos,_file.size());

    return newFile;
}

