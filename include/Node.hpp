#ifndef RNANUE_NODE_HPP
#define RNANUE_NODE_HPP

#include <iostream>
#include <vector>
#include <array>
#include <bitset>

class Interval {
    public:
        Interval(std::string chrom, char strand, std::string id, std::string name, int lower, int upper);
        ~Interval();
        // getter & setter
        int getLower() const;

        // operations
        std::string getAttribute(std::string key);



        void extend(int lower, int upper);

    private:
        std::string chrom;
        char strand;
        std::string id;
        std::string name;
        int lower;
        int upper;
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