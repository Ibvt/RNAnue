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
using RootNodes = std::vector<std::pair<std::string, Node*>>;

class IBPTree {
    public:
        IBPTree(po::variables_map params, int k);
        IBPTree();
        ~IBPTree();

        // constructs the IBPTree (either from annotations/clusters or both)
        void construct();
        void iterateFeatures(std::string featureFile);
        void insert(std::string chrom, const Interval& interval);

        std::map<std::string, std::string> getAttributes(const std::string& attributes);



        /*

        void iterateFeatures(std::string featureFile);
        void iterateClusters(std::string clusters);

        // operations on the tree structure
        void construct();
        void insert(std::string chrom, const Interval& interval);
        //void splitChild(InternalNode* parent, int index);
        void insertNonFull(Node* node, const Interval& interval);

        std::vector<Interval> search(const Interval& interval);
        void searchHelper(Node* node, const Interval& interval, std::vector<Interval>& result);
         */

    private:
        po::variables_map params;
        // tree structure for each chromosome
        RootNodes rootnodes; // generate a tree for each chromosome
        int order; // order of the IBPTree
};

#endif //RNANUE_IBTREE_HPP
