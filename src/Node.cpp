#include "Node.hpp"

Interval::Interval(int lower, int upper, char strand, std::string attributes) :
    lower(lower), upper(upper), strand(strand), attributes(attributes) {
}
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

