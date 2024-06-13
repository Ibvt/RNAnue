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
    iterateFeatures(params["features"].as<std::string>());
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
                                     std::make_pair(fields.start, fields.end), nullptr);
        } else {
            // only if the feature is an exon (needed for detect step)
            if(params["subcall"].as<std::string>() == "detect") {
                if(fields.type == "exon") {
                    // just make sure that its part of gene/transcript
                    if(intvl != nullptr) {
                        // just make sure that its part of the gene/transcript
                        if(intvl->isSubset(fields.start, fields.end)) {
                            std::string name = getTag(attr, std::vector<std::string>{"Parent"});
                            intvl->addJunction(name, std::make_pair(fields.start, fields.end));
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
    std::ifstream clusters(clusterFile);
    if(!clusters) {
        std::cout << helper::getTime() << " Cluster file " << clusterFile << " could not be opened!\n";
        EXIT_FAILURE;
    }

    IntervalData* firstSegment = nullptr;
    IntervalData* secondSegment = nullptr;

    // variables to buffer the current cluster
    std::string name = "";
    std::pair<std::string, std::string> chroms = std::make_pair("", "");
    std::pair<int, int> start = std::make_pair(0,0);
    std::pair<int, int> end = std::make_pair(0,0);
    std::pair<char, char> strand = std::make_pair(' ', ' ');
    std::pair<int, int> counts = std::make_pair(0,0);
    std::pair<int, int> size = std::make_pair(0,0);
    int nosplits = 0;

    IntervalData* intvl1 = nullptr;
    IntervalData* intvl2 = nullptr;

    std::string line;
    if(getline(clusters, line)) {} // ignore the first line (header)
    while(getline(clusters, line)) {
        std::string token;
        std::vector<std::string> tokens;
        std::istringstream ss(line);

        // split the input string
        while(getline(ss, token, '\t')) { tokens.push_back(token); }
        name = tokens[0];
        chroms = std::make_pair(tokens[1], tokens[5]);
        start = std::make_pair(std::stoi(tokens[3]), std::stoi(tokens[7]));
        end = std::make_pair(std::stoi(tokens[4]), std::stoi(tokens[8]));
        strand = std::make_pair(tokens[2][0], tokens[6][0]);
        nosplits = std::stoi(tokens[9]);

        if(nosplits < 2) { continue; } // ignore clusters with less than 2 splits reads
        std::vector<std::pair<Node*, IntervalData*>> fstOvlps = search(chroms.first,
                                                                       {start.first, end.first});

        if(fstOvlps.size() > 0) {
            for(auto& res : fstOvlps) {
                dtp::Interval tmpIntvl = res.second->getInterval();
                if(start.first < tmpIntvl.first && strand.first == res.second->getStrand()) {
                    res.second->setId(name + "/" + res.second->getId());
                    res.second->setName(name + "/" + res.second->getName());
                    res.second->setInterval({start.first, res.second->getInterval().second});
                    update(res.first->parent);
                }
                if(end.first > tmpIntvl.second && strand.first == res.second->getStrand()) {
                    res.second->setId(res.second->getId() + "/" + name);
                    res.second->setName(res.second->getName() + "/" + name);
                    res.second->setInterval({res.second->getInterval().first, end.first});
                    update(res.first->parent);
                }
                intvl1 = res.second;
            }
        } else {
            intvl1 = new IntervalData(chroms.second, strand.second, name, name, "cluster",
                                                   std::make_pair(start.second, end.second), nullptr);
            insert(*intvl1);
        }

        std::vector<std::pair<Node*, IntervalData*>> secOvlps = search(chroms.second,
                                                                       {start.second, end.second});
        if(secOvlps.size() > 0) {
            for(auto& res : secOvlps) {
                dtp::Interval tmpIntvl = res.second->getInterval();
                if(start.second < tmpIntvl.first && strand.second == res.second->getStrand()) {
                    res.second->setId(name + "/" + res.second->getId());
                    res.second->setName(name + "/" + res.second->getName());
                    res.second->setInterval({start.second, res.second->getInterval().second});
                    update(res.first->parent);
                }
                if(end.second > tmpIntvl.second && strand.second == res.second->getStrand()) {
                    res.second->setId(res.second->getId() + "/" + name);
                    res.second->setName(res.second->getName() + "/" + name);
                    res.second->setInterval({res.second->getInterval().first, end.second});
                    update(res.first->parent);
                }
                intvl2 = res.second;
            }
        } else {
            intvl2 = new IntervalData(chroms.second, strand.second, name, name, "cluster",
                                      std::make_pair(start.second, end.second), nullptr);
            insert(*intvl2);
        }
        intvl1->setSplit(intvl2);
        intvl2->setSplit(intvl1);

    }
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

std::vector<std::pair<Node*, IntervalData*>> IBPTree::search(std::string chrom, dtp::Interval interval) {
    Node* root = this->rootnodes[chrom]; // search for the root node
    std::vector<std::pair<Node*, IntervalData*>> result;
    searchIter(root, interval, result);

    std::sort(result.begin(), result.end());
    auto last = std::unique(result.begin(), result.end());
    result.erase(last, result.end());
    return result;
}

void IBPTree::searchIter(Node* node, const dtp::Interval& interval,
                         std::vector<std::pair<Node*, IntervalData*>>& results) {
    if(node->isLeaf) {
        for(int i=0; i<node->keys.size();++i) {
            if(isOverlapping(interval, node->keys[i].first)) {
                results.push_back(std::make_pair(node, node->keys[i].second));
            }
        }

        if(node->next != nullptr) {
            if(isOverlapping(interval, node->next->keys[0].first)) {
                searchIter(node->next, interval, results);
            }
        }
    } else {
        int i = 0;
        while(i < node->keys.size() &&
                (interval.first > node->keys[i].first.first &&
                        (interval.first > node->keys[i].first.second))){
            i++;
        }
        if(node->children[i] != nullptr) {
            searchIter(node->children[i], interval, results);
        }
    }
}

void IBPTree::update(Node* node) {
    // update parent of node (until root is reached)
    if(node == nullptr) {
        return;
    } else {
        for(int i=0; i<node->keys.size(); ++i) {
            dtp::Interval intvl = node->keys[i].first;
            std::cout << node->children[i]->keys.size() << "\n";
            for(int j=0; j<node->children[i]->keys.size(); ++j) {
                if(node->children[i]->keys[j].first < intvl) {
                    intvl.first = node->children[i]->keys[j].first.first;
                }
                if(node->children[i]->keys[j].first.second > intvl.second) {
                    intvl.second = node->children[i]->keys[j].first.second;
                }
            }
            node->keys[i].first = intvl;
        }
        update(node->parent);
    }
}

bool IBPTree::isOverlapping(dtp::Interval intvl1, dtp::Interval intvl2) {
    dtp::Interval intvl = {std::max(intvl1.first, intvl2.first),
                           std::min(intvl1.second, intvl2.second)};
    if(intvl.first <= intvl.second) {
        return true;
    } else {
        return false;
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
