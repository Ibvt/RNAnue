#include "Node.hpp"

Interval::Interval(int lower, int upper, char strand, std::string attributes) :
    lower(lower), upper(upper), strand(strand), attributes(attributes) {
}
Interval::~Interval() {}

void Interval::extend(int lower, int upper) {
    this->lower = lower;
    this->upper = upper;
}

// create new node
Node::Node(int order) {
    this->keys = new int[order]; // initialize keys array
    this->children = new Node*[order+1]; // array of pointers to children
}

