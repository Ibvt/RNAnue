#include "IBPTree.hpp"

IBPTree::IBPTree(po::variables_map params, int k) : params(params) {
    this->rootnodes = std::vector<std::pair<std::string, Node*>>();
    this->order = k; // can be later extracted from config file

    std::cout << helper::getTime() << " Constructing IBPTree\n";
    construct();
}

IBPTree::IBPTree() {}
IBPTree::~IBPTree() {}


// constructs the IBPTree (either from annotations/clusters or both)
void IBPTree::construct() {
    std::string subcall = params["subcall"].as<std::string>();
    // fill tree with annotation
    if(subcall == "detect") { // in the case of detect only annotation is needed
        iterateFeatures(params["features"].as<std::string>());
    }
}

void IBPTree::iterateFeatures(std::string featuresFile) {
    // iterate over features
    std::ifstream gff(featuresFile);
    if(!gff) {
        std::cout << helper::getTime() << " Annotation file " << featuresFile << " could not be opened!\n";
        EXIT_FAILURE;
    }

    // null pointer
    Interval* intvl = nullptr;

    // create variables for GFF file
    std::string chrom = "", source = "", biotype = "", start = "";
    std::string end = "", score = "", attributes = "";
    char strand = ' ';

    std::string line;
    while(getline(gff, line)) {
        if(line[0] == '#') { // ignore header lines
            continue;
        }

        std::string token;
        std::vector<std::string> tokens;
        std::istringstream ss(line);

        // split the input string
        while(getline(ss, token, '\t')) {
            tokens.push_back(token);
        }
        chrom = tokens[0]; // chromosome
        source = tokens[1]; // source
        biotype = tokens[2]; // biotype
        start = std::stoi(tokens[3]); // start
        end = std::stoi(tokens[4]); // end
        score = tokens[5]; // score
        strand = tokens[6].toChar();
        attributes = tokens[8]; // attributes

        if(biotype == "region") {
            continue;
        }
        std::map<std::string, std::string> attr = getAttributes(attributes);
        // check if id is in the attributes - can also be substring

        if(attr.find("ID") != attr.end()) {
            std::string id = attr["ID"];
        } else {
            std::string id = "";
        }

        if(intvl == nullptr) {
            std::string id = "", name = "";
            // the next element that is read should be 'gene'
            if(attr.find("gene_name") != attr.end()) {
                name = attr["gene_name"];
            } else {
                if(attr.find("Name") != attr.end()) {
                    std::string name = attr["Name"];
                }
            }
            intvl = new Interval(chrom, strand, id, name, start, end);
        }
    }
    gff.close();
}

void IBPTree::insert(std::string chrom, const Interval& interval) {


}

// return
std::map<std::string, std::string> IBPTree::getAttributes(const std::string& attributes) {
    std::map<std::string, std::string> attr; // output to be map with key-value pairs
    std::istringstream ss(attributes); // create string stream
    std::string token;
    while(getline(ss, token, ';')) {
        std::string key = token.substr(0, token.find("="));
        std::string value = token.substr(token.find("=") + 1);
        attr.insert(std::make_pair(key, value));
    }
    return attr;
}





/*
void IBPTree::insert(std::string chrom, const Interval& interval) {
    // check if chrom is already in the tree (and if not add it)
    auto it = std::find_if(rootnodes.begin(), rootnodes.end(),
                           [&chrom](const std::pair<std::string, Node*>& p) {
        return p.first == chrom;
    });
    if(it == rootnodes.end()) {
        rootnodes.push_back(std::make_pair(chrom, nullptr));
    } else {
        if(it->second == nullptr) { // the tree is empty
            it->second = new Node(order);





        }
    }










        // search for the node and insert the interval
        if(it->second == nullptr) {
            it->second = new LeafNode(interval);
        } else {
            if(it->second->isLeaf()) {
                LeafNode* leaf = static_cast<LeafNode*>(it->second);
                if(leaf->intervals.size() == 2 * order - 1) {
                    InternalNode* newRoot = new InternalNode();
                    newRoot->children.push_back(it->second);
                    splitChild(newRoot, 0);
                    it->second = newRoot;
                }
            }
        }
    }







    if(root == nullptr) { // tree is empty
        root = new LeafNode(interval);
        return;
    } else{
        if(root->isLeaf()) {
            LeafNode* found = search(chrom);
            LeafNode* leaf = static_cast<LeafNode*>(root);
            if(leaf->intervals.size() == 2 * order - 1) {
                InternalNode* newRoot = new InternalNode();
                newRoot->children.push_back(root);
                splitChild(newRoot, 0);
                root = newRoot;
            }
        }
    }
    // make sure the root is a leaf node
    if (root->isLeaf()) {
        LeafNode* leaf = static_cast<LeafNode*>(root);
        if (leaf->intervals.size() == 2 * order - 1) {
            InternalNode* newRoot = new InternalNode();
            newRoot->children.push_back(root);
            splitChild(newRoot, 0);
            root = newRoot;
        }
    }
    insertNonFull(root, interval);
}

// Function to split a child of a node
void IBPTree::splitChild(InternalNode* parent, int index) {
    InternalNode* child = static_cast<InternalNode*>(parent->children[index]);
    InternalNode* newChild = new InternalNode();

    parent->keys.insert(parent->keys.begin() + index, child->keys[degree - 1]);
    parent->children.insert(parent->children.begin() + index + 1, newChild);

    for (int j = 0; j < degree - 1; j++) {
        newChild->keys.push_back(child->keys[j + degree]);
    }
    child->keys.resize(degree - 1);

    for (int j = 0; j < degree; j++) {
        newChild->children.push_back(child->children[j + degree]);
    }
    child->children.resize(degree);
}

void IBPTree::insertNonFull(Node* node, const Interval& interval) {
    if (node->isLeaf()) {
        LeafNode* leaf = static_cast<LeafNode*>(node);
        leaf->keys.push_back(interval);
        std::sort(leaf->keys.begin(), leaf->keys.end(), [](const Interval& a, const Interval& b) {
            return a.low < b.low;
        });
    } else {
        InternalNode* internal = static_cast<InternalNode*>(node);
        int i = internal->keys.size() - 1;
        while (i >= 0 && interval.low < internal->keys[i]) {
            i--;
        }
        i++;
        if (internal->children[i]->isLeaf() && static_cast<LeafNode*>(internal->children[i])->keys.size() == 2 * degree - 1) {
            splitChild(internal, i);
            if (interval.low > internal->keys[i])
                i++;
        }
        insertNonFull(internal->children[i], interval);
    }
}

// Function to search for intervals overlapping with the given interval
std::vector<Interval> IBPTree::search(const Interval& interval) {
    std::vector<Interval> result;
    if (root != nullptr)
        searchHelper(root, interval, result);
    return result;
}

// Helper function for searching intervals
void IBPTree::searchHelper(Node* node, const Interval& interval, std::vector<Interval>& result) {
    if (node->isLeaf()) {
        LeafNode* leaf = dynamic_cast<LeafNode*>(node);
        for (const Interval& key : leaf->keys) {
            if (interval.low <= key.high && interval.high >= key.low)
                result.push_back(key);
        }
    } else {
        InternalNode* internal = dynamic_cast<InternalNode*>(node);
        int i = 0;
        while (i < internal->keys.size() && interval.low > internal->keys[i])
            i++;
        searchHelper(internal->children[i], interval, result);
    }
}

void IBPTree::iterateClusters(std::string clusterFile) {
    // iterate over clusters
    std::ifstream clusters(clusterFile);
    if(!clusters) {
        std::cout << helper::getTime() << " Cluster file " << clusters << " could not be opened!\n";
        EXIT_FAILURE;
    }

    std::pair<int, int> start = std::make_pair(0, 0);
    std::pair<int, int> end = std::make_pair(0, 0);
    std::pair<char, char> strand = std::make_pair(' ', ' ');
    std::pair<int, int> counts = std::make_pair(0, 0);
    std::pair<int, int> size = std::make_pair(0, 0);

    std::string line;
    while(getline(clusters, line)) {
        std::string token;
        std::vector<std::string> tokens;
        std::istringstream ss(line);

        // split the input string
        while(getline(ss, token, '\t')) {
            tokens.push_back(token);
        }

        start = std::make_pair(std::stoi(tokens[1]), std::stoi(tokens[2]));
        end = std::make_pair(std::stoi(tokens[3]), std::stoi(tokens[4]));
        strand = std::make_pair(tokens[5][0], tokens[6][0]);
        counts = std::make_pair(std::stoi(tokens[7]), std::stoi(tokens[8]));
        size = std::make_pair(std::stoi(tokens[9]), std::stoi(tokens[10]));

        insert(tokens[0], Interval(start, end, strand, counts, size));
    }
}


*/