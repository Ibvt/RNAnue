//
// Created by Richard Albin Schaefer on 1/25/24.
//

#ifndef SEQRICKSHAW_HPP
#define SEQRICKSHAW_HPP

#include <iostream>
#include <tuple>
#include <cstddef>
#include <algorithm>
#include <future>
#include <mutex>
#include <thread>
#include <chrono>
#include <map>

#include<omp.h>

// boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

// seqan
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>

#include "StateTransition.hpp"
#include "Preproc.hpp"
#include "Utility.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

class SeqRickshaw {
    public:
        SeqRickshaw(po::variables_map params);
        ~SeqRickshaw();

        void start(pt::ptree sample);
        void distribute(dtp::SubPathsMap reads); // process the reads
        void processing(dtp::PathVector subreads);

        // the workers that distribute the processing of the reads
        void workerSE(seqan3::sequence_file_input<> &fin, seqan3::sequence_file_output<> &fout, std::mutex outMutex);

        // functions involve in the splitting procedures
        int getLinesPerThread(fs::path path);
        dtp::SubPathsMap split(fs::path forward, fs::path reverse);
        void splitcall(fs::path inFile, fs::path outDir, int linesPerThread);

        void testfun();

        // output

    private:
        po::variables_map params;
        StateTransition adpt5Tables;
        StateTransition adpt3Tables;

        int modetrm;
        int threads;
};

#endif //SEQRICKSHAW_HPP
