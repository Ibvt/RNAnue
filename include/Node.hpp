#ifndef RNANUE_NODE_HPP
#define RNANUE_NODE_HPP

#include <iostream>
#include <vector>

class Interval {
    public:
        Interval();
        ~Interval();

        void extend();

    private:
        int lower;
        int upper;
        char strand;
        std::string attributes;
};

class Node {
    public:
        virtual ~Node() {}
        virtual bool isLeaf() const = 0;
};

class InternalNode: public Node {
    public:
        InternalNode();
        ~InternalNode();

        bool isLeaf();

        std::vector<int> keys;
        std::vector<Node*> children;
};

class LeafNode: public Node {
    public:
        LeafNode(const Interval& interval);
        ~LeafNode();

        bool isLeaf() const;

    private:
        std::vector<Interval> intervals;
        std::string attributes;
};

#endif //RNANUE_NODE_HPP