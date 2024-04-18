#ifndef RNANUE_NODE_HPP
#define RNANUE_NODE_HPP

#include <iostream>
#include <vector>
#include <array>
#include <bitset>

class Interval {
    public:
        Interval();
        Interval(int lower, int upper, char strand, std::string attributes);
        ~Interval();



        // getter & setter
        int getLower() const;


        void extend(int lower, int upper);

    private:
        int lower;
        int upper;
        char strand;
        std::string attributes;
};

class Node {
    public:
        Node(int order);
        void addInterval(Interval& interval);

    private:
        int order;
        std::vector<int> keys;
        std::vector<Node*> children;
        std::bitset<1> leaf;
};

#endif //RNANUE_NODE_HPP