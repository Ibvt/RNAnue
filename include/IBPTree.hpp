#ifndef RNANUE_IBTREE_HPP
#define RNANUE_IBTREE_HPP

// Standard
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>

// Boost
#include <boost/program_options.hpp>

// Class
#include "Utility.hpp"
#include "Node.hpp"
#include "DataTypes.hpp"

namespace po = boost::program_options;
using RootNodes = std::map<std::string, Node*>;

class IBPTree {
    public:
        IBPTree(po::variables_map params, int k);
        IBPTree();
        ~IBPTree();

        // constructs the IBPTree (either from annotations/clusters or both)
        void construct();
        void iterateFeatures(std::string featureFile);
        void iterateClusters(std::string clusterFile);

        // tree operations
        Node* getRoot(IntervalData& data);
        void insert(IntervalData& data);
        void insertIter(Node* node, IntervalData& data);
        void splitNode(Node* node, int index);
        std::vector<IntervalData*> search(std::string chrom, dtp::Interval interval);
        void searchIter(Node* node, const dtp::Interval& interval, std::vector<IntervalData*> results);
        bool isOverlapping(dtp::Interval intvl1, dtp::Interval intvl2);


        std::map<std::string, std::string> getAttributes(std::string& attributes);
        std::string getTag(std::map<std::string, std::string> attributes, const std::vector<std::string>& keys);

        // add tree operations
        void printTree();
        void traverse(Node* parent, Node* child, int link, std::ostream& oss) const;

    private:
        po::variables_map params;
        // tree structure for each chromosome
        RootNodes rootnodes; // generate a tree for each chromosome
        int order; // order of the IBPTree
};

#endif //RNANUE_IBTREE_HPP
