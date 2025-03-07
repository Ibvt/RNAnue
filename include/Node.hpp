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
                 std::string biotype, std::string product, dtp::Interval interval,
                 IntervalData* split);
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
        std::string getProduct() const;
        void setProduct(std::string product);
        dtp::Interval getInterval();
        void setInterval(dtp::Interval interval);
        dtp::SpliceJunctions getJunctions() const;
        void setJunctions(dtp::SpliceJunctions junctions);
        IntervalData* getSplit() const;
        void setSplit(IntervalData* split);

        // operations
        void addJunction(std::string name, std::pair<size_t,size_t> junction);
        bool isSubset(int start, int end);
        bool isOverlapping(dtp::Interval intvl1, dtp::Interval intvl2);
        void printNode();

    private:
        std::string chrom;
        char strand;
        std::string id;
        std::string name;
        std::string biotype;
        std::string product;
        dtp::Interval interval;
        dtp::SpliceJunctions junctions;
        IntervalData* split;
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