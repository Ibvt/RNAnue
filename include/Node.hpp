#ifndef RNANUE_NODE_HPP
#define RNANUE_NODE_HPP

#include <iostream>
#include <vector>
#include <array>

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
        Node(int order);

    private:
        int* keys;
        Node** children;
};

#endif //RNANUE_NODE_HPP