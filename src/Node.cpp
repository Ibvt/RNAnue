#include "Node.hpp"

IntervalData::IntervalData(std::string chrom, char strand, std::string id, std::string name,
                   std::string biotype, dtp::Interval interval) :
    chrom(chrom), strand(strand), id(id), name(name), biotype(biotype), interval(interval) {
}
IntervalData::~IntervalData() {}
bool IntervalData::operator>(const dtp::Interval& other) const { return this->interval.first > other.first; }
bool IntervalData::operator<(const dtp::Interval& other) const { return this->interval.first < other.first; }

// getter & setter
std::string IntervalData::getChrom() const { return chrom; }
void IntervalData::setChrom(std::string chrom) { this->chrom = chrom; }
char IntervalData::getStrand() const { return strand; }
void IntervalData::setStrand(char strand) { this->strand = strand; }
std::string IntervalData::getName() const { return name; }
void IntervalData::setName(std::string name) { this->name = name; }
std::string IntervalData::getBiotype() const { return biotype; }
void IntervalData::setBiotype(std::string biotype) { this->biotype = biotype; }
dtp::Interval IntervalData::getInterval() { return interval; }
void IntervalData::setInterval(dtp::Interval interval) { this->interval = interval; }
std::string IntervalData::getJunctions() const { return junctions; }
void IntervalData::setJunctions(std::string junctions) { this->junctions = junctions; }

// operations
void IntervalData::addJunction(std::string junction) { this->junctions += junction; }
//bool IntervalData::isSubset(int start, int end) { return (start >= this->lower && end >= this->upper); }
void IntervalData::printNode() {
    std::cout << "----------\nChrom: " << this->chrom << "\nStrand: " << this->strand;
    std::cout << "\nBnds: [" << this->interval.first << "," << this->interval.second << "]";
    std::cout << "\nID: " << this->id << "\nName: " << this->name << "\nBiotype: " << this->biotype;
    std::cout << "\nJunctions: " << this->junctions << "\n----------";
}
bool IntervalData::isSubset(int start, int end) {
    return (start >= this->interval.first && end <= this->interval.second);
}


// create new node
Node::Node(int k) : order{k}, keys{}, children{}, next{nullptr}, parent{nullptr}, isLeaf{false} {}

// operations
void Node::addInterval(IntervalData& data) {
    // add Interval to the node (by comparing the interval) - should be sorted operator > in IntervalData
    int i = 0;
    while(i < keys.size() && data > keys[i].first) { i++; }
    this->keys.insert(this->keys.begin() + i, {data.getInterval(), &data});
}

void Node::addKey(dtp::Interval interval, int index) {
    this->keys.insert(this->keys.begin() + index, {interval, nullptr});
}

dtp::Interval Node::calcNewKey() {
    dtp::Interval intvl = {std::string::npos, 0};
    for(int i = 0; i < keys.size(); i++) {
        if(keys[i].first.first < intvl.first) { intvl.first = keys[i].first.first; }
        if(keys[i].first.second > intvl.second) { intvl.second = keys[i].first.second; }
    }
    return intvl;
}

void Node::addChild(Node* child) {
    this->children.push_back(child);
}

std::string Node::keysToString() {
    std::string out = "";
    for(int i = 0; i < keys.size(); i++) {
        if(i == 0) { out += "|";}
        out += "[" + std::to_string(keys[i].first.first) + "," + std::to_string(keys[i].first.second) + "]|";
    }
    return out;
}

Node* Node::getChild(int index) {
    return this->children[index];
}
