#include "IBPTree.hpp"

IBPTree::IBPTree(po::variables_map params, int k) : params(params) {
    this->rootnodes = std::map<std::string, Node*>();
    this->order = k; // can be later extracted from config file

    std::cout << helper::getTime() << "Constructing IBPTree(s) using ";
    std::cout << "features from " << params["features"].as<std::string>() << "\n";
    construct();
    std::cout << helper::getTime() << "IBPTree(s) constructed...\n";

    printTree(); // save tree to file
    std::cout << helper::getTime() << "... and saved to file...\n";
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
        std::cout << helper::getTime() << "Annotation file " << featuresFile << " could not be opened!\n";
        EXIT_FAILURE;
    }

    // null pointer
    IntervalData* intvl = nullptr;
    dtp::FeatureFields fields; // stores the current fields of the features
    std::string line;
    std::string junction = "";

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
        fields.seqid = tokens[0];
        fields.source = tokens[1];
        fields.type = tokens[2];
        fields.start = std::stoi(tokens[3]);
        fields.end = std::stoi(tokens[4]);
        fields.score = tokens[5];
        fields.strand = tokens[6][0];
        fields.phase = tokens[7][0];
        fields.attributes = tokens[8];

        // ignore regions (as this would be return as overlap with all intervals)
        if(fields.type == "region") {
            continue;
        }
        std::map<std::string, std::string> attr = getAttributes(fields.attributes); // get attributes
        // get tags needed for the intervals
        std::string id = getTag(attr, std::vector<std::string>{"ID"});
        std::string name = getTag(attr, std::vector<std::string>{"gene_name", "Name"});
        std::string biotype = getTag(attr, std::vector<std::string>{"gene_biotype", "gene_type"});

        if(attr.find("Parent") == attr.end()) { // new gene/transcript indicated by missing 'Parent' attribute
            if(intvl != nullptr) { // there exists and interval object (add to tree)
                insert(*intvl);
                intvl = nullptr;
            }
            intvl = new IntervalData(fields.seqid, fields.strand, id, name, biotype,
                                     std::make_pair(fields.start, fields.end));
        } else {
            // only if the feature is an exon (needed for detect step)
            if(params["subcall"].as<std::string>() == "detect") {
                if(fields.type == "exon") {
                    // just make sure that its part of gene/transcript
                    if(intvl != nullptr) {
                        // just make sure that its part of the gene/transcript
                        if(intvl->isSubset(fields.start, fields.end)) {
                            junction = std::to_string(fields.start) + "-" + std::to_string(fields.end);
                            intvl->addJunction(junction);
                        }
                    }
                }
            }
        }
    }
    // add the last interval
    if(intvl != nullptr) {
        insert(*intvl);
    }
    gff.close();
}

void IBPTree::iterateClusters(std::string clusterFile) {
    /*
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
    }*/
}

Node* IBPTree::getRoot(IntervalData& data) {
    Node* root;
    if(this->rootnodes.find(data.getChrom()) == this->rootnodes.end()) { // chrom-tree is emptry (create new rootnode)
        root = new Node(this->order);
        root->isLeaf = true;
        this->rootnodes.insert(std::make_pair(data.getChrom(), root));
    } else {
        root = this->rootnodes[data.getChrom()];
    }
    return root;
}

void IBPTree::insert(IntervalData &data) {
    Node* root = getRoot(data);
    insertIter(root, data); // insert the new interval
    if(root->keys.size() == this->order) {
        Node* newRoot = new Node(this->order);
        newRoot->children.push_back(root);
        splitNode(newRoot, 0);
        root = newRoot;
        this->rootnodes[data.getChrom()] = root;
    }
}

void IBPTree::insertIter(Node* node, IntervalData& data) {
    if(node->isLeaf) {
        node->addInterval(data);
    } else {
        int childnum = 0;
        while(childnum < node->keys.size() && data > node->keys[childnum].first) {
            childnum++;
        }
        insertIter(node->children[childnum], data);

        if(node->children[childnum]->keys.size() == order) {
            splitNode(node, childnum);
        }
    }
}

// split the node
void IBPTree::splitNode(Node* parent, int index) {
   // std::cout << "Splitting node " << parent->keysToString() << " at index " << index << "\n";
    Node* child = parent->getChild(index);
    Node* newChild = new Node(this->order);
    int mid = ((this->order+2-1)/2); // value for order=6

    newChild->isLeaf = child->isLeaf;
    newChild->keys.assign(child->keys.begin() + mid, child->keys.end());
    child->keys.resize(this->order-1);

    // update parent
    parent->children.insert(parent->children.begin() + index + 1, newChild);
    parent->keys.insert(parent->keys.begin() + index, {child->calcNewKey(), nullptr});

    if(child->isLeaf) {
        newChild->next = child->next;
        child->next = newChild;
    } else {
        newChild->children.assign(child->children.begin() + mid, child->children.end());
        child->children.resize(mid + 1);
    }
}

// get attribute from the attributes fields
std::string IBPTree::getTag(std::map<std::string, std::string> attributes, const std::vector<std::string>& keys) {
    std::string element = "";
    for(auto& key : keys) {
        if(attributes.find(key) != attributes.end()) {
            return attributes[key];
        }
    }
    return element;
}

std::map<std::string, std::string> IBPTree::getAttributes(std::string& attributes) {
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

// tree operations
void IBPTree::printTree() {
    fs::path outDir = params["outdir"].as<std::string>();
    fs::path treeOutDir = outDir / fs::path("IBPTree");

    for(auto& root : this->rootnodes) {
        fs::path treeOutFile = treeOutDir / fs::path(root.first + ".txt");
        std::ofstream treeOutStream(treeOutFile.string(), std::ios::trunc);
        if (treeOutStream.is_open()) {
            traverse(nullptr, root.second, -1, treeOutStream);
            treeOutStream.close();
        } else {
            std::cerr << helper::getTime() << "Error: Unable to open the file " << treeOutFile.string() << ".\n";
        }
    }
}

void IBPTree::traverse(Node* parent, Node* child, int link, std::ostream& ofs) const {
    if(parent != nullptr && child != nullptr) {
        ofs << parent->keysToString() << "\tlink(" << link << ")\t" << child->keysToString() << "\n";
    }
    for(int i=0; i<child->children.size(); i++) {
        traverse(child, child->children[i], i, ofs);
    }
}

/*
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



*/
/*
if(root->getLeaf()) { // root is a leaf node
    root->addInterval(data);
    if(root->getKeys().size() > this->order-1) {
        splitNode(root);
    }
} else {
    dtp::Interval interval = std::make_pair(5300,5400);
    Node* child = traverse(interval, root);
    if(child == nullptr) {
        std::cout << "nullptr" << std::endl;
    } else {
        std::cout << "new Child: " << child->keysToString() << std::endl;
    }

*/




/*
std::cout << "Child: " << child->keysToString() << std::endl;
 */

