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

// Class
#include "Utility.hpp"

namespace po = boost::program_options;
namespace pt = boost::property_tree;
namespace fs = boost::filesystem;

struct Segment {
    std::string refid;
    uint32_t flag;
    uint32_t start;
    uint32_t end;

    Segment() : refid(""), flag(0), start(0), end(0) {};
    Segment(std::string _refid, uint32_t _flag, unsigned int _start, unsigned int _end) :
            refid(_refid), flag(_flag), start(_start), end(_end) {};
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
    Clustering();

    void iterate(std::string splits);
    void overlaps(std::vector<Cluster> &clusterlist);
    bool startPosCmp(Cluster &a, Cluster &b);
    void start(pt::ptree sample);
    void sumup();
};




#endif //RNANUE_CLUSTERING_HPP
