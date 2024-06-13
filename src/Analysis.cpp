#include "Analysis.hpp"

//
Analysis::Analysis(po::variables_map _params) : params(_params), repcount{0}, readcount{0} {
    int k = 3; // can be later extracted from config file
    // extract features (e.g., IBPTree)
    this->features = IBPTree(params, k); // create IBPT w/ annotations
    if(params["clust"].as<std::bitset<1>>() == std::bitset<1>("1")) {
        std::cout << helper::getTime() << "Add clusters to IBPTree\n";
        fs::path outDir = params["outdir"].as<std::string>();
        fs::path clusterFile = outDir / fs::path("clustering") / fs::path("clusters.tab");
        this->features.iterateClusters(clusterFile.string());
    }
    this->freq = std::map<std::pair<std::string, std::string>, double>();
    this->suppreads = std::map<dtp::IntKey, std::vector<double>>();
    this->complementarities = std::map<dtp::IntKey, std::vector<std::vector<double>>>();
    this->hybenergies = std::map<dtp::IntKey, std::vector<std::vector<double>>>();
}

Analysis::~Analysis() {
}

void Analysis::start(pt::ptree sample, pt::ptree condition) {
    this->freq.clear(); // reset freq
    // retrieve input and output files
    pt::ptree input = sample.get_child("input");
    std::string single = input.get<std::string>("single");
    std::string splits = input.get<std::string>("splits");

    pt::ptree output = sample.get_child("output");
    std::string interactions = output.get<std::string>("interactions");

    //
    this->freq = std::map<std::pair<std::string, std::string>, double>();
    this->condition = condition.get<std::string>("condition"); // store current condition

    this->readcount = 0;

    processSingle(single);
    processSplits(splits, interactions);

    this->repcount++; // increase replicate count

    normalize();
}

void Analysis::processSingle(std::string& single) {
    seqan3::sam_file_input inputSingle{single}; // match reads
    std::vector<size_t> ref_lengths{}; // header
    for(auto &info : inputSingle.header().ref_id_info) {
        ref_lengths.push_back(std::get<0>(info));
    }
    std::deque<std::string> refIds = inputSingle.header().ref_ids();
    std::optional<int32_t> refOffset = std::nullopt;
    std::string orient = ""; // e.g., sense (ss) or antisense (as)
    char strand = ' ';
    uint32_t start = 0, end = 0;
    for(auto& rec : inputSingle) {
        refOffset = rec.reference_position();
        start = rec.reference_position().value();
        end = rec.reference_position().value() + rec.sequence().size();
        if(static_cast<bool>(rec.flag() & seqan3::sam_flag::on_reverse_strand)) { strand = '-'; } else { strand = '+'; }

        std::vector<std::pair<Node*, IntervalData*>> ovlps = features.search(refIds[rec.reference_id().value()],
                                                                             {start, end});
        if(ovlps.size() > 0) {
            for(auto& ovl : ovlps) {
                if(strand == ovl.second->getStrand()) { orient = "sense"; } else { orient = "antisense"; }
                addToFreqMap(std::make_pair(orient, ovl.second->getName()), 1.0/(double)ovlps.size()); // add to freq
            }
        }
        readcount++;
    }
}

void Analysis::processSplits(std::string& splits, std::string& output) {
    seqan3::sam_file_input inputSplits{splits}; // match reads
    std::ofstream outInts(output, std::ios::out); // output interactions
    writeInteractionsHeader(outInts); // write header in output file

    std::vector<size_t> ref_lengths{}; // header
    for(auto &info : inputSplits.header().ref_id_info) {
        ref_lengths.push_back(std::get<0>(info));
    }
    std::deque<std::string> refIds = inputSplits.header().ref_ids();

    // variables
    std::string outLine = ""; // buffers the output file for interaction entry
    std::bitset<1> skip = 0; // skip if no match
    char strand = ' '; // strand of segment
    std::string orient = ""; // orientation of segment (ss or as)
    uint32_t start = 0, end = 0; // start and end of segment
    std::string refName = ""; // reference name

    // information from the tags
    seqan3::sam_tag_dictionary tags;
    std::pair<float, float> filters; // complementarity and hybridization energy
    std::vector<std::string> aligns; // alignment of complementarity

    // segment information (e.g., gene name, product, etc.) - keys are the as or ss
    std::map<std::string, std::vector<IntervalData*>> segmentInfo = {};

    for(auto&& rec : inputSplits | seqan3::views::chunk(2)) {
        outLine = (*rec.begin()).id(); // QNAME
        // reset
        skip = 0; // reset skip
        segmentInfo.clear(); // reset segment info
        filters = std::make_pair(0.0, 0.0); // reset filters
        aligns.clear();

        // retrieve complementarity and energy
        tags = (*rec.begin()).tags();
        filters = std::make_pair(std::get<float>(tags["XC"_tag]), std::get<float>(tags["XE"_tag]));

        dtp::IntKey intKey = {}; // key to store the interaction (partner1, orientation1, partner2, orientation2)

        for (auto &split: rec) {
            if (static_cast<bool>(split.flag() & seqan3::sam_flag::on_reverse_strand)) {
                strand = '-';
            } else { strand = '+'; }
            outLine += "\t" + std::string(1, strand); // strand
            start = split.reference_position().value();
            end = split.reference_position().value() + split.sequence().size();
            outLine += "\t" + std::to_string(start) + "\t" + std::to_string(end); // start & end
            refName = refIds[split.reference_id().value()];
            outLine += "\t" + refName; // reference name

            tags = split.tags(); // retrieve the alignment information
            aligns.push_back(std::get<std::string>(tags["XA"_tag]));

            // match with features
            std::vector<std::pair<Node *, IntervalData *>> ovlps = features.search(refIds[split.reference_id().value()],
                                                                                   {start, end});
            segmentInfo = {}; // reset segmentInfo
            if(ovlps.size() > 0) {
                segmentInfo.insert({"sense", {}});
                segmentInfo.insert({"antisense", {}});
                for (auto &ovl: ovlps) {
                    if (strand == ovl.second->getStrand()) { orient = "sense"; } else { orient = "antisense"; }
                    segmentInfo[orient].push_back(ovl.second); // add to segmentInfo
                    addToFreqMap(std::make_pair(orient, ovl.second->getName()), 1.0/(double)ovlps.size()); // add to freq


                    /*
                    std::string info = ovl.second->getName() + "\t" + ovl.second->getBiotype() + "\t" +
                                       ovl.second->getStrand() + "\t" + "." + "\t" + orient;
                    segmentInfo[orient].push_back(info);*/
                }
                // check which one to select
                if (segmentInfo["sense"].size() == 1) {
                    IntervalData* data = segmentInfo["sense"][0];
                    outLine += "\t" + data->getName() + "\t" + data->getBiotype() + "\t" + data->getStrand();
                    outLine += "\t.\tsense";
                    intKey.push_back(dtp::IntPartner(data->getName(),"sense"));
                } else {
                    if (segmentInfo["antisense"].size() == 1) {
                        IntervalData* data = segmentInfo["antisense"][0];
                        outLine += "\t" + data->getName() + "\t" + data->getBiotype() + "\t" + data->getStrand();
                        outLine += "\t.\tsense";
                        intKey.push_back(dtp::IntPartner(data->getName(),"antisense"));
                    } else {
                        skip = 1; // skip if not unique
                    }
                }
            } else { skip = 1; } // skip if no match
        }
        outLine += "\t" + std::to_string(filters.first); // complementarity
        outLine += "\t" + aligns[0] + "\t" + aligns[1]; // alignment
        outLine += "\t" + std::to_string(filters.second) + "\n"; // energy
        if (skip == 0) {
            outInts << outLine;
            addToSuppReads(intKey[0], intKey[1]);
            addToCompl(intKey[0], intKey[1], filters.first);
            addToHybEnergies(intKey[0], intKey[1], filters.second);
        } // write to file (if not skipped)
        readcount++;
    }
    outInts.close();
}

void Analysis::writeInteractionsHeader(std::ofstream& fout) {
    fout << "qname\tfst_seg_strd\tfst_seg_strt\tfst_seg_end\tfst_seg_ref\tfst_seg_name\tfirst_seg_bt\t";
    fout << "fst_seg_anno_strd\tfst_seg_prod\tfst_seg_ori\tsec_seg_strd\tsec_seg_strt\tsec_seg_end\t";
    fout << "sec_seg_ref\tsec_seg_name\tsec_seg_bt\tsec_seg_anno_strd\tsec_seg_prod\tsec_seg_ori\t";
    fout << "cmpl\tfst_seg_compl_aln\tsec_seg_cmpl_aln\tmfe\n";
}

void Analysis::writeAllIntsHeader(std::ofstream& fout) { fout << "fst_rna\tsec_rna\tfst_rna_ori\tsec_rna_ori";}
void Analysis::addToAllIntsHeader(std::ofstream &fout, std::string key) { fout << "\t" << key; }

void Analysis::normalize() {
    // sum up all frequencies
    double sum = 0.0;
    for(auto& f : this->freq) {
        sum += f.second;
    }

    // normalize
    for(auto& f : this->freq) {
        f.second = f.second/sum;
    }

    // print
    /*
    for(auto& f : freq) {
        std::cout << f.first.first << " " << f.first.second << " " << f.second << std::endl;
    }*/
}

void Analysis::addToFreqMap(std::pair<std::string, std::string> key, double value) {
if(freq.count(key) == 0) {
        freq.insert(std::make_pair(key, value));
    } else {
        freq[key] += value;
    }
}

bool operator<(const dtp::IntKey& lhs, const dtp::IntKey& rhs) {
    return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

void Analysis::addToSuppReads(dtp::IntPartner p1, dtp::IntPartner p2) {
    dtp::IntKey key;
    if(p1 < p2) { key = {p1, p2}; } else { key = {p2, p1}; }
    // check if key in map
    if(suppreads.find(key) == suppreads.end()) {
        std::vector<double> suppReadsRepl = {}; // new list for all replicates (up to the current one)
        for(int i=0;i<repcount;++i) { suppReadsRepl.push_back(0.0); }
        suppReadsRepl.push_back(1.0); // add 1.0 for the last (current) replicate
        suppreads.insert({key, suppReadsRepl});
    } else { // key found in map
        std::vector<double>& suppReadsRepl = this->suppreads[key];
        if(suppReadsRepl.size() < repcount) { // not yet at the current replicate (probably skipped a few)
            for(int i=suppReadsRepl.size(); i<repcount; ++i) {
                suppReadsRepl.push_back(0.0);
            }
            suppReadsRepl.push_back(1.0);
        } else {
            if(suppReadsRepl.size() > repcount) {
                suppReadsRepl[repcount] += 1.0;
            } else {
                suppReadsRepl.push_back(1.0);
            }
        }
    }
}

void Analysis::addToCompl(dtp::IntPartner p1, dtp::IntPartner p2, double complementarity) {
    dtp::IntKey key;
    if(p1 < p2) { key = {p1, p2}; } else { key = {p2, p1}; }
    // check if key in map
    if(complementarities.find(key) == complementarities.end()) {
        std::vector<std::vector<double>> CmplsRepl = {}; // new list for all replicates (up to the current one)
        for(int i=0;i<repcount;++i) { CmplsRepl.push_back({}); }
        CmplsRepl.push_back({complementarity}); // add 1.0 for the last (current) replicate
        complementarities.insert({key, CmplsRepl});
    } else {
        std::vector<std::vector<double>>& cmplsRepl = this->complementarities[key];
        if(cmplsRepl.size() < repcount) { // not yet at the current replicate (probably skipped a few)
            for(int i=cmplsRepl.size(); i<repcount; ++i) { cmplsRepl.push_back({}); }
            cmplsRepl.push_back({complementarity});
        } else {
            if(cmplsRepl.size() > repcount) {
                cmplsRepl[repcount].push_back(complementarity);
            } else {
                cmplsRepl.push_back({complementarity});
            }
        }
    }
}

void Analysis::addToHybEnergies(dtp::IntPartner p1, dtp::IntPartner p2, double hybenergy) {
    dtp::IntKey key;
    if(p1 < p2) { key = {p1, p2}; } else { key = {p2, p1}; }
    // check if key in map
    if(hybenergies.find(key) == hybenergies.end()) {
        std::vector<std::vector<double>> hybRepl = {}; // new list for all replicates (up to the current one)
        for(int i=0;i<repcount;++i) { hybRepl.push_back({}); }
        hybRepl.push_back({hybenergy}); // add 1.0 for the last (current) replicate
        hybenergies.insert({key, hybRepl});
    } else {
        std::vector<std::vector<double>>& hybRepl = this->hybenergies[key];
        if(hybRepl.size() < repcount) { // not yet at the current replicate (probably skipped a few)
            for(int i=hybRepl.size(); i<repcount; ++i) { hybRepl.push_back({}); }
            hybRepl.push_back({hybenergy});
        } else {
            if(hybRepl.size() > repcount) {
                hybRepl[repcount].push_back(hybenergy);
            } else {
                hybRepl.push_back({hybenergy});
            }
        }
    }
}

double Analysis::calcGCS(std::vector<double>& complementarities) {
    if(complementarities.size() == 0) { return 0.0; }
    double median = stats::median(complementarities);
    double maximum = *std::max_element(complementarities.begin(), complementarities.end());
    return median * maximum;
}

double Analysis::calcGHS(std::vector<double>& hybenergies) {
    if(hybenergies.size() == 0) { return 1000.0; }
    double median = stats::median(hybenergies);
    double minimum = *std::min_element(hybenergies.begin(), hybenergies.end());
    if(median < 0 && minimum < 0) {
        return -sqrt(median * minimum);
    } else {
        if(median > 0 && minimum > 0) {
            return sqrt(median * minimum);
        } else {
            return sqrt(abs(median * minimum));
        }
    }
}

double Analysis::calcStat(dtp::IntKey key, int x) {
    std::pair<std::string, std::string> p1 = {key[0].orientation, key[0].partner};
    std::pair<std::string, std::string> p2 = {key[1].orientation, key[1].partner};
    double p;
    if(freq.find(p1) == freq.end() || freq.find(p2) == freq.end()) {
        p = 0.0;
    } else {
        if(p1 == p2) {
            p = freq[p1]*freq[p2];
        } else {
            p = 2*freq[p1]*freq[p2];
        }
    }
    ma::binomial_distribution<double> binomial(this->readcount, p);
    double pval = 1.0 - ma::cdf(binomial, x);
    return pval;
}

void Analysis::writeAllInts() {
    std::string outDirPath = params["outdir"].as<std::string>();
    fs::path outDirDetectPath = fs::path(outDirPath) / fs::path("analysis");
    fs::path allIntsPath = fs::path(outDirDetectPath) / fs::path("allints.txt");
    std::ofstream fout(allIntsPath.string(), std::ios::out);

    if(fout.is_open()) {
        for(auto& entry : suppreads) {
            fout << entry.first[0].partner << "\t" << entry.first[1].partner;
            fout << "\t" << entry.first[0].orientation << "\t" << entry.first[1].orientation;
            for(int i=0;i<entry.second.size();++i) {
                fout << "\t" << entry.second[i]; // counts
                double gcs = calcGCS(complementarities[entry.first][i]); // complementarity score (GCS)
                if(gcs != 0.0) { fout << "\t" << gcs; } else { fout << "\t."; }
                double ghs = calcGHS(hybenergies[entry.first][i]); // Hybridization score (GHS)
                if(ghs != 1000.0) { fout << "\t" << ghs; } else { fout << "\t."; }
                // statistical
                // call with key(RNA1/2 and orient) and counts
                double stat = calcStat(entry.first, entry.second[i]);
                fout << "\t" << stat;
            }
            if(entry.second.size() < repcount) {
                for(int i=entry.second.size();i<repcount;++i) {
                    fout << "\t0\t.\t.\t.";
                }
            }
            fout << "\n";
        }
        fout.close();
        std::cout << helper::getTime() << "Interactions written to file: " << allIntsPath.string() << "\n";
    }
}

