// Standard
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// Boost
#define BOOST_TEST_MAIN
#define BOOST_THROW_EXCEPTIONS
#include <boost/test/unit_test.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// Class
#include "Base.hpp"
#include "Data.hpp"

// namespace
namespace po = boost::program_options;
namespace fs = boost::filesystem;

bool compareFiles(const fs::path& p1, const fs::path& p2) {
    std::ifstream ifs1(p1.string());
    std::ifstream ifs2(p2.string());

    std::istream_iterator<char> b1(ifs1), e1;
    std::istream_iterator<char> b2(ifs2), e2;

    BOOST_CHECK_EQUAL_COLLECTIONS(b1, e1, b2, e2);
}

BOOST_AUTO_TEST_CASE(DataHandlingPreprocSE) {
    fs::path testdir{fs::path(__FILE__).parent_path()};
    fs::path outdir{testdir / "data" / "datahandling" / "results_SE"}; // test-outdir
    fs::path ctrls{testdir / "data" / "datahandling" / "rawreads_SE" / "ctrls"}; // test-ctrls
    fs::path trtms{testdir / "data" / "datahandling" / "rawreads_SE" / "trtms"}; // test-trtms

    // Inserting dummy key-value pairs
    po::variables_map params;
    params.insert(std::make_pair("ctrls", po::variable_value(ctrls.string(), false)));
    params.insert(std::make_pair("trtms", po::variable_value(trtms.string(), false)));
    params.insert(std::make_pair("outdir", po::variable_value(outdir.string(), false)));
    params.insert(std::make_pair("readtype", po::variable_value(std::string("SE"), false)));
    params.insert(std::make_pair("subcall", po::variable_value(std::string("preproc"), false)));

    // Create Data Object
    Data data(params);

    // compare files
    fs::path ref{testdir / "data" / "datahandling" / "dataPreprocSE.json"}; // ground truth
    fs::path res{outdir / "data.json"}; // generated data file
    compareFiles(ref, res);
}

BOOST_AUTO_TEST_CASE(DataHandlingPreprocPE) {
    fs::path testdir{fs::path(__FILE__).parent_path()};
    fs::path outdir{testdir / "data" / "datahandling" / "results_PE"}; // test-outdir
    fs::path ctrls{testdir / "data" / "datahandling" / "rawreads_PE" / "ctrls"}; // test-ctrls
    fs::path trtms{testdir / "data" / "datahandling" / "rawreads_PE" / "trtms"}; // test-trtms

    // Inserting dummy key-value pairs
    po::variables_map params;
    params.insert(std::make_pair("ctrls", po::variable_value(ctrls.string(), false)));
    params.insert(std::make_pair("trtms", po::variable_value(trtms.string(), false)));
    params.insert(std::make_pair("outdir", po::variable_value(outdir.string(), false)));
    params.insert(std::make_pair("readtype", po::variable_value(std::string("PE"), false)));
    params.insert(std::make_pair("subcall", po::variable_value(std::string("preproc"), false)));

    // Create Data Object
    Data data(params);

    // compare files
    fs::path ref{testdir / "data" / "datahandling" / "dataPreprocPE.json"}; // ground truth
    fs::path res{outdir / "data.json"}; // generated data file
    compareFiles(ref, res);
}
