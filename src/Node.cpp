#include "Node.hpp"

Interval::Interval(int lower, int upper, char strand, std::string attributes) :
    lower(lower), upper(upper), strand(strand), attributes(attributes) {
}
Interval::~Interval() {}

void Interval::extend(int lower, int upper) {
    this->lower = lower;
    this->upper = upper;
}

InternalNode::InternalNode() {
}

InternalNode::~InternalNode() {
}

LeafNode::LeafNode(const Interval& interval) {
    intervals.push_back(interval);
}

bool LeafNode::isLeaf() const {
    return true;
}

LeafNode::~LeafNode() {
}





