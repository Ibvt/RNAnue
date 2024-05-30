#ifndef RNANUE_NODE_HPP
#define RNANUE_NODE_HPP

// Standard
#include <iostream>
#include <vector>
#include <array>
#include <bitset>

// Class
#include "DataTypes.hpp"

class IntervalData {
    public:
        IntervalData(std::string chrom, char strand, std::string id, std::string name,
                 std::string biotype, dtp::Interval interval);
        ~IntervalData();

        // operator overloading
        bool operator>(const dtp::Interval& other) const;
        bool operator<(const dtp::Interval& other) const;

        // getter & setter
        std::string getChrom() const;
        void setChrom(std::string chrom);
        char getStrand() const;
        void setStrand(char strand);
        std::string getId() const;
        void setId(std::string id);
        std::string getName() const;
        void setName(std::string name);
        std::string getBiotype() const;
        void setBiotype(std::string biotype);
        dtp::Interval getInterval();
        void setInterval(dtp::Interval interval);
        std::string getJunctions() const;
        void setJunctions(std::string junctions);

        // operations
        void addJunction(std::string junction);
        bool isSubset(int start, int end);
        void printNode();

    private:
        std::string chrom;
        char strand;
        std::string id;
        std::string name;
        std::string biotype;
        dtp::Interval interval;
        std::string junctions;
};

class Node {
    public:
        // class variables
        int order;
        std::vector<std::pair<dtp::Interval, IntervalData*>> keys;
        std::vector<Node*> children;
        Node* next; // link to the next node
        Node* parent; // link to the parent node
        bool isLeaf;

        // constructor
        Node(int order);

        // operations
        void addInterval(IntervalData& data);
        void addKey(dtp::Interval interval, int index);
        dtp::Interval calcNewKey();
        void addChild(Node* child);
        Node* getChild(int index);
        std::string keysToString();

};

#endif //RNANUE_NODE_HPP