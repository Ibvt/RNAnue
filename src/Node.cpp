#include "Node.hpp"

Interval::Interval(std::string chrom, char strand, std::string id, std::string name, int lower, int upper) :
    chrom(chrom), strand(strand), id(""), name(""), lower(-1), upper(-1)  {}

Interval::~Interval() {}

void Interval::extend(int lower, int upper) {
    this->lower = lower;
    this->upper = upper;
}

int Interval::getLower() const {
    return lower;
}

// create new node
Node::Node(int order) : order(order), keys(), children() {}

void Node::addInterval(Interval& interval) {
    // check if the node is a leaf
    if(leaf == std::bitset<1>("1")) {
        // check if there is space in the node
        if(keys.size() < order) {

        }

    }

    // check if there is still space in the node
    if(keys.size() < order) { // node is not full


    }
}

