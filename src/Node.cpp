#include "Node.hpp"

Interval::Interval(std::string chrom, char strand, std::string id, std::string name,
                   std::string biotype, int lower, int upper) :
    chrom(chrom), strand(strand), id(""), name(""), biotype(biotype), lower(-1), upper(-1)  {}

Interval::~Interval() {}

void Interval::narrow(int lower, int upper) {
    this->lower = lower;
    this->upper = upper;
}

// getter & setter
int Interval::getLower() const {
    return lower;
}

void Interval::setId(std::string id) {
    this->id = id;
}
std::string Interval::getId() const {
    return id;
}

std::string Interval::getChrom() const {
    return chrom;
}

void Interval::setJunction(std::pair<int,int> junction) {
    junctions.push_back(junction);
}


bool Interval::isSubset(int start, int end) {
    return (start >= this->lower && end >= this->upper);
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

