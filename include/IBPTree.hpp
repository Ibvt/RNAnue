#ifndef RNANUE_IBTREE_HPP
#define RNANUE_IBTREE_HPP

// Standard
#include <iostream>
#include <fstream>
#include <vector>

// Boost
#include <boost/program_options.hpp>

// Class
#include "Utility.hpp"
#include "Node.hpp"

namespace po = boost::program_options;

class IBPTree {
    public:
        IBPTree(po::variables_map params, int order);
        ~IBPTree();

        void iterateFeatures(std::string featureFile);
        void iterateClusters(std::string clusters);

        // operations on the tree structure
        void construct();
        void insert(std::string chrom, const Interval& interval);
        void splitChild(InternalNode* parent, int index);
        void insertNonFull(Node* node, const Interval& interval);

        std::vector<Interval> search(const Interval& interval);
        void searchHelper(Node* node, const Interval& interval, std::vector<Interval>& result);

    private:
        po::variables_map params;
        // tree structure for each chromosome
        std::vector<std::pair<std::string,Node*>> rootnodes; // gener
        int order; // order of the IBPTree
};

#endif //RNANUE_IBTREE_HPP
