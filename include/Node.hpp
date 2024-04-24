#ifndef RNANUE_NODE_HPP
#define RNANUE_NODE_HPP

#include <iostream>
#include <vector>
#include <array>
#include <bitset>

class Interval {
    public:
        Interval(std::string chrom, char strand, std::string id, std::string name,
                 std::string biotype, int lower, int upper);
        ~Interval();

        // getter & setter
        int getLower() const;
        void setId(std::string id);
        std::string getId() const;
        std::string getChrom() const;
        void setJunction(std::pair<int,int> junction);

        // operations
        std::string getAttribute(std::string key);
        bool isSubset(int start, int end);
        void narrow(int lower, int upper);

    private:
        std::string chrom;
        char strand;
        std::string id;
        std::string name;
        std::string biotype;
        int lower;
        int upper;
        std::vector<std::pair<int,int>> junctions;
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