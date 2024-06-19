#include "Data.hpp"

Data::Data(po::variables_map params) : params(params) {
    std::string subcall = params["subcall"].as<std::string>();
    fs::path outDir = fs::path(params["outdir"].as<std::string>());

    // create output directory (e.g., outdir)
    helper::createOutDir(outDir, std::cout);
    // create directory for subcall (e.g., outdir/preproc)
    fs::path outSubcallDir = outDir / fs::path(subcall);
    helper::createDir(outSubcallDir, std::cout);
    helper::createDir(outDir / fs::path("tmp"), std::cout);

    // preprocessing
    if(subcall == "preproc") { preprocDataPrep(); }
    if(subcall == "align") { alignDataPrep(); }
    if(subcall == "detect") {
        fs::path outIBTree = outDir / fs::path("IBPTree");
        helper::createDir(outIBTree, std::cout);
        detectDataPrep();
    }
    if(subcall == "clustering") { clusteringDataPrep();}
    if(subcall == "analysis") { analysisDataPrep(); }
}

Data::~Data() {
    fs::path outDir = fs::path(params["outdir"].as<std::string>());
    fs::path tmpDir = outDir / fs::path("tmp");
    helper::deleteDir(tmpDir);
}

void Data::preprocDataPrep() {
    // retrieve paths that contain the reads
    fs::path ctrlsPath = fs::path(this->params["ctrls"].as<std::string>());
    fs::path trtmsPath = fs::path(this->params["trtms"].as<std::string>());

    GroupsPath groups = getGroupsPath(ctrlsPath, trtmsPath);
    getCondition(groups);
}

void Data::alignDataPrep() {
    if(params["preproc"].as<std::bitset<1>>() == std::bitset<1>("0")) {
        std::cout << helper::getTime() << "Skipping preprocessing\n";
        preprocDataPrep(); // same data collecting as in preprocessing
    } else {
        // make sure that data has been preprocessed (or at least selected)
        if(params["preproc"].as<std::bitset<1>>() == std::bitset<1>("1")) {
            fs::path ctrlsPath = fs::path(params["outdir"].as<std::string>()) / "preproc/ctrls";
            fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / "preproc/trtms";

            /*
            * trtms -> "/Users/..../trtms"
            * (ctrls -> "/Users/.../crtls")
            */
            GroupsPath groups = getGroupsPath(ctrlsPath, trtmsPath);
            getCondition(groups);
        }
    }
}

void Data::detectDataPrep() {
    fs::path ctrlsPath = fs::path(params["outdir"].as<std::string>()) / "align/ctrls";
    fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / "align/trtms";

    GroupsPath groups = getGroupsPath(ctrlsPath, trtmsPath);
    getCondition(groups);
}

void Data::clusteringDataPrep() {
fs::path ctrlsPath = fs::path(params["outdir"].as<std::string>()) / "detect/ctrls";
    fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / "detect/trtms";

    GroupsPath groups = getGroupsPath(ctrlsPath, trtmsPath);
    getCondition(groups);
}

void Data::analysisDataPrep() {
    fs::path ctrlsPath = fs::path(params["outdir"].as<std::string>()) / "detect/ctrls";
    fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / "detect/trtms";

    GroupsPath groups = getGroupsPath(ctrlsPath, trtmsPath);
    getCondition(groups);
}

GroupsPath Data::getGroupsPath(fs::path& ctrls, fs::path& trtms) {
    GroupsPath groups;
    std::cout << helper::getTime() << "Retrieve the path that contain the reads\n";

    if(trtms != "") { // '--trtms' has been set in the config file or cmdline
        if(ctrls != "") { // '--ctrls' has been set in the config or cmdline
            groups.insert(std::make_pair("ctrls", ctrls));
        } else { // '--ctrls' has not been set in either config or cmdline - still continue
            std::cout << helper::getTime() << "### WARNING - '--ctrls' has not been set. ";
            std::cout << "This runs RNAnue without control data\n";
        }
        groups.insert(std::make_pair("trtms", trtms));
    } else { // abort - trtms have not been set in either config file or cmdline
        std::cerr << "### ERROR - RNAnue is aborted! No treatment data specified" << std::endl;
        exit(EXIT_FAILURE);
    }
    return groups;
}

//
void Data::getCondition(GroupsPath& groups) {
    // ptree object containing
    pt::ptree ptData, ptSubcall, ptPath, ptGroup, ptCondition;

    for(auto& group : groups) { // iterate through the groups (e.g., ctrls, trtms)
        std::cout << helper::getTime() << group.first << ": " << group.second << " ";
        if(fs::exists(group.second)) { // check if the provided path exists...
            std::cout << "has been found in the filesystem!\n";
            if(fs::is_directory(group.second)) { // .. and check if its a directory (and not a file)
                // retrieve content (conditions) within group (e.g., rpl-exp,...) as absolute path
                PathVector conditions = helper::listDirFiles(group.second);

                for(auto& condition : conditions) {
                    if(fs::exists(condition) && !fs::is_directory(condition)) {
                        std::cout << helper::getTime() << "### WARNING - " << condition << " is not a directory! ";
                        std::cout << "This condition will be skipped!\n";
                        continue;
                    }

                    ptCondition = getData(group.first, condition);
                    ptGroup.push_back(std::make_pair("", ptCondition));

                }
                ptSubcall.add_child(group.first, ptGroup);
                ptGroup.erase("");

            } else {
                std::cout << helper::getTime() << "### ERROR - " << group.second << " is not a directory! ";
                exit(EXIT_FAILURE);
            }
        } else {
            std::cout << helper::getTime() << "### ERROR - " << group.second << " has not been found in the filesystem!\n";
            exit(EXIT_FAILURE);
        }
    }

    // write data structure to file
    std::string subcall = params["subcall"].as<std::string>();
    dataStructure.add_child(subcall, ptSubcall);
    fs::path dataStructurePath = fs::path(params["outdir"].as<std::string>()) / fs::path("data.json");
    jp::write_json(dataStructurePath.string(), dataStructure);
    std::cout << helper::getTime() << "Data structure has been written to " << dataStructurePath << "\n";
}

//
pt::ptree Data::getData(std::string group, fs::path& condition) {
    std::string subcall = params["subcall"].as<std::string>();
    std::string outdir = params["outdir"].as<std::string>();

    pt::ptree ptCondition, ptFiles, ptSamples, ptSample, ptOutput;

    // scan the condition directory for files
    PathVector dataFiles = helper::listDirFiles(condition);
    dataFiles = helper::filterDirFiles(dataFiles, subcall);

    /*
    for(auto& file : dataFiles) {
        std::cout << file << std::endl;
    }*/

    // define helper variables
    int numberElements = getNumberElements(dataFiles);
    std::vector<std::string> sampleKeys = getSampleKeys();

    fs::path conditionOutDir = fs::path(outdir);
    conditionOutDir /= fs::path(subcall);
    conditionOutDir /= fs::path(group);
    conditionOutDir /= fs::path(condition.filename());

    int elementCounter = 0;
    for(unsigned i=0;i<dataFiles.size();++i) {
        ptFiles.put(sampleKeys[i % numberElements], dataFiles[i].string());
        if(elementCounter == numberElements-1) {
            ptSample.add_child("input", ptFiles);
            // determine output TODO
            ptOutput = getOutputData(ptSample, conditionOutDir);
            ptSample.add_child("output", ptOutput);

            ptSamples.push_back(std::make_pair("", ptSample));
            ptSample.erase("input");
            ptSample.erase("output");
            for(unsigned j=0;j<sampleKeys.size();++j) {
                ptFiles.erase(sampleKeys[j]);
            }
            elementCounter = 0;
        } else {
            ++elementCounter;
        }
    }
    ptCondition.put("condition", condition.stem().string());
    ptCondition.add_child("samples", ptSamples);

    return ptCondition;
}


int Data::getNumberElements(PathVector& vec) {
    std::string subcall = params["subcall"].as<std::string>();
    int numberElements;

    if(subcall == "preproc") { numberElements = (params["readtype"].as<std::string>() == "PE") ? 2 : 1; }
    if(subcall == "align") { numberElements = (params["readtype"].as<std::string>() == "PE") ? 5 : 1; }
    if(subcall == "detect") {
        std::vector<std::string> keys = {"preproc_matched", "R1only_matched", "R2only_matched", "unmerged_matched"};
        numberElements = 1;
        // for now only consider preproc_matched
        // TODO: additional files when unmerged/unfiltered reads are used (paired-end)
    }
    if(subcall == "clustering") { numberElements = 3; }
    if(subcall == "analysis") { numberElements = 3; }
    return numberElements;
}

std::vector<std::string> Data::getSampleKeys() {
    std::string subcall = params["subcall"].as<std::string>();
    std::vector<std::string> sampleKeys;
    if(subcall == "preproc") { sampleKeys = {"forward", "reverse"}; }

    if(subcall == "align") {
        if(params["readtype"].as<std::string>() == "SE") {
            sampleKeys = {"forward"};
        } else {
            if(params["readtype"].as<std::string>() == "PE") {
                sampleKeys = {"forward", "R1only", "R1unmerged", "R2only", "R2unmerged"};
            }
        }
    }

    if(subcall == "detect") { sampleKeys = {"matched"}; }
    if(subcall == "clustering") { sampleKeys = {"multsplits", "single", "splits"}; }
    if(subcall == "analysis") { sampleKeys = {"multsplits", "single", "splits"}; }
    return sampleKeys;
}

pt::ptree Data::getOutputData(pt::ptree& input, fs::path& conditionOutDir) {
    std::string subcall = params["subcall"].as<std::string>();

    pt::ptree output;
    if(subcall == "preproc") { output = getPreprocOutputData(input, conditionOutDir); }
    if(subcall == "align") { output = getAlignOutputData(input, conditionOutDir); }
    if(subcall == "detect") { output = getDetectOutputData(input, conditionOutDir); }
    if(subcall == "analysis") { output = getAnalysisOutputData(input, conditionOutDir); }
    return output;
}

pt::ptree Data::getPreprocOutputData(pt::ptree& input, fs::path& conditionOutDir) {
    pt::ptree output;

    // replace input path with output path (results/...)
    fs::path inForward = fs::path(input.get<std::string>("input.forward"));
    std::string outForward = helper::replacePath(conditionOutDir, inForward).string();

    // create output files using the input files with an additonal suffix
    std::string outPreproc = helper::addSuffix(outForward,"_preproc", {"_R1","fwd"});
    output.put("forward", outPreproc); // write output back to ptree

    // using paired-end reads results in additional files for unassembled reads
    if(params["readtype"].as<std::string>() == "PE") {
        // replace
        fs::path inReverse = fs::path(input.get<std::string>("input.reverse"));
        std::string outReverse = helper::replacePath(conditionOutDir, inReverse).string();

        // create outfiles using the input files with an additional suffix
        std::string r1only = helper::addSuffix(outForward,"_R1only", {"_R1","fwd"});
        std::string r2only = helper::addSuffix(outReverse,"_R2only", {"_R2","rev"});

        std::string r1unmerged = helper::addSuffix(outForward, "_R1unmerged", {"_R1","fwd"});
        std::string r2unmerged = helper::addSuffix(outReverse, "_R2unmerged", {"_R2","rev"});

        output.put("R1only", r1only); // push both back to output ptree
        output.put("R2only", r2only);

        output.put("R1unmerged", r1unmerged);
        output.put("R2unmerged", r2unmerged);
    }

    return output;
}

pt::ptree Data::getAlignOutputData(pt::ptree& input, fs::path& conditionOutDir) {
    pt::ptree output;

    // replace input path with output path (results/...)
    fs::path inPreproc = fs::path(input.get<std::string>("input.forward"));
    std::string outPreproc = helper::replacePath(conditionOutDir,
                                              inPreproc.replace_extension(".bam")).string();

    // create output (before adding suffix for file)
    output.put("matched", helper::addSuffix(outPreproc, "_matched", {}));

    if(params["readtype"].as<std::string>() == "PE") {
        // create output file R1 only reads (
        fs::path inR1only = fs::path(input.get<std::string>("input.R1only"));
        std::string outR1only = helper::replacePath(conditionOutDir,
                                                    inR1only.replace_extension(".bam")).string();
        output.put("matched_R1only", helper::addSuffix(outR1only, "_matched", {}));

        // create output file R2 only reads
        fs::path inR2only = fs::path(input.get<std::string>("input.R2only"));
        std::string outR2only = helper::replacePath(conditionOutDir,
                                                    inR1only.replace_extension(".bam")).string();
        output.put("matched_R2only", helper::addSuffix(outR2only, "_matched", {}));

        // create out file for unmerged reads
        output.put("matched_unmerged", helper::addSuffix(outPreproc, "_unmerged_matched",
                                                 {"_preproc"}));
    }

    return output;
}

pt::ptree Data::getDetectOutputData(pt::ptree& input, fs::path& conditionOutDir) {
    pt::ptree output;

    // replace input path with output path (results/...)
    fs::path inMatched = fs::path(input.get<std::string>("input.matched"));
    std::string outMatched = helper::replacePath(conditionOutDir, inMatched).string();
    output.put("single", helper::addSuffix(outMatched, "_single", {"_matched"}));
    output.put("splits", helper::addSuffix(outMatched, "_splits", {"_matched"}));
    output.put("multsplits", helper::addSuffix(outMatched, "_multsplits", {"_matched"}));
    return output;
}

pt::ptree Data::getAnalysisOutputData(pt::ptree& input, fs::path& conditionOutDir) {
    pt::ptree output;

    // replace input path with output path (results/...)
    fs::path inSplits = fs::path(input.get<std::string>("input.splits"));
    std::string outSplits = helper::replacePath(conditionOutDir, inSplits.replace_extension(".txt")).string();
    output.put("interactions", helper::addSuffix(outSplits, "_interactions", {"_splits"}));
    return output;
}

template <typename Callable>
void Data::callInAndOut(Callable f) {

    // retrieve paths for parameters
    fs::path outDir = fs::path(params["outdir"].as<std::string>());
    std::string subcallStr = params["subcall"].as<std::string>();
    fs::path outSubcallDir = outDir / fs::path(subcallStr);

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
        if(params["subcall"].as<std::string>() != "clustering") {
            // create directory for groups (e.g., ctrls, trtms)
            helper::createDir(outGroupDir, std::cout); // not needed for clustering
        }

        conditions = subcall.get_child(groups[i]);
        BOOST_FOREACH(pt::ptree::value_type const &v, conditions.get_child("")) {
            pt::ptree condition = v.second;

            // create directory for condition (e.g., rpl_exp)
            if(params["subcall"].as<std::string>() != "clustering") {
                outConditionDir = outGroupDir / fs::path(condition.get<std::string>("condition"));
                helper::createDir(outConditionDir, std::cout);
            }

            samples = condition.get_child("samples");
            // iterate over samples
            BOOST_FOREACH(pt::ptree::value_type const &w, samples.get_child("")) {
                pt::ptree sample = w.second;
                f(sample, condition);
            }
        }
    }
}

//
void Data::preproc() {
    std::cout << helper::getTime() << "Start the Preprocessing\n";
    SeqRickshaw srs(params);
    callInAndOut(std::bind(&SeqRickshaw::start, srs, std::placeholders::_1));
}

void Data::align() {
    std::cout << helper::getTime() << "Start the Alignment\n";
    Align aln(params);
    callInAndOut(std::bind(&Align::start, aln, std::placeholders::_1));
}

void Data::detect() {
    std::cout << helper::getTime() << "Start the Split Read Calling\n";
    SplitReadCalling src(params);
    callInAndOut(std::bind(&SplitReadCalling::start, &src, std::placeholders::_1, std::placeholders::_2));
    src.writeStats();
}

void Data::clustering() {
    std::cout << helper::getTime() << "Start the Clustering of the Split Reads\n";
    Clustering clu(params);
    callInAndOut(std::bind(&Clustering::start, &clu, std::placeholders::_1));
    clu.sumup();
}

void Data::analysis() {
    std::cout << helper::getTime() << "Start the Analysis of the Split Reads\n";
    Analysis ana(params);
    callInAndOut(std::bind(&Analysis::start, &ana, std::placeholders::_1, std::placeholders::_2));
    // write file output files
    ana.writeAllInts();
    ana.writeAllIntsCounts();
    ana.writeAllIntsJGF();
    ana.writeStats();
}

