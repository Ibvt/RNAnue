#ifndef RNANUE_CLUSTERING_HPP
#define RNANUE_CLUSTERING_HPP

// Standard
#include <algorithm>
#include <omp.h>

// Boost
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>

// SeqAn3
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/utility/views/chunk.hpp>

// Class
#include "Utility.hpp"

namespace po = boost::program_options;
namespace pt = boost::property_tree;
namespace fs = boost::filesystem;

struct Segment {
    std::string refid;
    seqan3::sam_flag flag;
    uint32_t start;
    uint32_t end;

    Segment() : refid(""), flag(seqan3::sam_flag{}), start(0), end(0) {};
    Segment(std::string refid, seqan3::sam_flag flag, unsigned int start, unsigned int end) :
            refid(refid), flag(flag), start(start), end(end) {};
};

struct Cluster {
    std::vector<Segment> elements;
    int count;

    Cluster() : elements({}), count(1) {};
    bool operator<(const Cluster& a) const {
        return elements[0].start < a.elements[0].start;
    }
};


class Clustering {
private:
    po::variables_map params;
    std::vector<Cluster> result;

public:
    Clustering(po::variables_map params);

    void iterate(std::string splits);
    void overlaps(std::vector<Cluster> &clusterlist);
    bool startPosCmp(Cluster &a, Cluster &b);
    void start(pt::ptree sample);
    void sumup();
};




#endif //RNANUE_CLUSTERING_HPP
