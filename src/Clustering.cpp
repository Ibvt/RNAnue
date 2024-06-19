#include <iostream>
#include "Clustering.hpp"

Clustering::Clustering(po::variables_map params) : params(params), result({}) {}

void Clustering::iterate(std::string splits) {
    seqan3::sam_file_input fin{splits};
    int chunks = 5;
    std::vector<Cluster> subset;
    // header
    std::vector<size_t> ref_lengths{};
    for(auto &info : fin.header().ref_id_info) {
        ref_lengths.push_back(std::get<0>(info));
    }
    std::deque<std::string> ref_ids = fin.header().ref_ids();
    for(auto&& rec: fin | seqan3::views::chunk(2)) {
        Cluster cl;
        for(auto& split : rec) { // push
            //seqan3::debug_stream << split << std::endl;
            std::optional<int32_t> refId = split.reference_id();
            seqan3::sam_flag flag = split.flag(); // SAMFLAG
            uint32_t start = split.reference_position().value();
            uint32_t end = start + split.sequence().size()-1;
            Segment seg(ref_ids[refId.value()], flag, start, end);
            cl.elements.push_back(seg);
        }

        // always have the left segment occurring first
        if(cl.elements[0].start > cl.elements[1].start) {
            Segment cl0 = cl.elements[0];
            cl.elements[0] = cl.elements[1];
            cl.elements[1] = cl0;
        }
        subset.push_back(cl);
        if(subset.size() == chunks) {
            overlaps(subset);
            result.insert(result.end(),subset.begin(),subset.end());
            subset.clear();
        }
    }
    if(subset.size() > 0) {
        overlaps(subset);
        result.insert(result.end(),subset.begin(),subset.end());
        //  std::cout << "results size: " << result.size() << std::endl;
        subset.clear();
    }
    overlaps(result);
    // print clusters
    /*
    std::cout << "clusters\n";
    for(auto& cl : result) {
        std::cout << cl.elements[0].refid << "\t" << cl.elements[0].start << "\t" << cl.elements[0].end << "\t";
        std::cout << cl.elements[1].refid << "\t" << cl.elements[1].start << "\t" << cl.elements[1].end << "\t";
        std::cout << cl.count << std::endl;
    }*/
}

//
bool Clustering::startPosCmp(Cluster &a, Cluster &b) {
    return a.elements[0].start < b.elements[0].start;
}

void Clustering::sumup() {
    // retrieve output directory
    fs::path output = fs::path(params["outdir"].as<std::string>()) / fs::path("clustering");
    fs::path clusterResults = output / fs::path("clusters.tab");
    std::ofstream outputFile(clusterResults.string());

    std::cout << helper::getTime() << "Write clusters to file: " << clusterResults << std::endl;

    outputFile << "clustID\tfst_seg_chr\tfst_seg_strd\tfst_seg_strt\tfst_seg_end\t"
                  "sec_seg_chr\tsec_seg_strd\tsec_seg_strt\tsec_seg_end\tno_splits\t"
                  "fst_seg_len\tsec_seg_len\n";

    int clustID = 1;
    for(unsigned i=0;i<result.size();++i) {
        if(outputFile.is_open()) {
            outputFile << "cluster" << clustID++ << "\t";
            outputFile << result[i].elements[0].refid << "\t";
            if(!static_cast<bool>(result[i].elements[0].flag & seqan3::sam_flag::on_reverse_strand)) {
                outputFile << "+" << "\t";
            } else {
                outputFile << "-" << "\t";
            }
            outputFile << result[i].elements[0].start << "\t";
            outputFile << result[i].elements[0].end << "\t";

            outputFile << result[i].elements[1].refid << "\t";
            if(!static_cast<bool>(result[i].elements[0].flag & seqan3::sam_flag::on_reverse_strand)) {
                outputFile << "+" << "\t";
            } else {
                outputFile << "-" << "\t";
            }
            outputFile << result[i].elements[1].start << "\t";
            outputFile << result[i].elements[1].end << "\t";
            outputFile << result[i].count << "\t";
            outputFile << (result[i].elements[0].end+1) - result[i].elements[0].start << "\t";
            outputFile << (result[i].elements[1].end+1) - result[i].elements[1].start << "\n";
        }
    }
    outputFile.close();

    /*
    for(unsigned i=0;clusters.size();++i) {
        //std::cout << clusters[i].count << std::endl;
    }*/
}
void Clustering::overlaps(std::vector<Cluster> &clusterlist) {
    // sort by first segment (and the second)
    std::sort(clusterlist.begin(), clusterlist.end(), [this](Cluster &a, Cluster &b) {
        if(a.elements[0].start == b.elements[0].start) {
            return a.elements[1].start < b.elements[1].start;
        }
        return a.elements[0].start < b.elements[0].start;
    });

    int clustdist = params["clustdist"].as<int>();

    uint32_t s1Start, s1End, s2Start, s2End;
    uint32_t xs1Start, xs1End, xs2Start, xs2End;

    for(unsigned i=0;i<clusterlist.size();++i) {
        for(unsigned j=i+1; j<clusterlist.size();++j) {
            // check if the next cluster is too far away
            if(clusterlist[j].elements[0].start > clusterlist[i].elements[0].end + clustdist) {
                break;
            }

            uint32_t s1Start = clusterlist[i].elements[0].start+1;
            uint32_t s1End = clusterlist[i].elements[0].end+1;
            uint32_t s2Start = clusterlist[i].elements[1].start+1;
            uint32_t s2End = clusterlist[i].elements[1].end+1;

            uint32_t xs1Start = clusterlist[j].elements[0].start+1;
            uint32_t xs1End = clusterlist[j].elements[0].end+1;
            uint32_t xs2Start = clusterlist[j].elements[1].start+1;
            uint32_t xs2End = clusterlist[j].elements[1].end+1;

            // checks if the first overlaps
            if((xs1Start <= s1End + clustdist) && (xs1Start >= s1Start - clustdist)) {
                if((s2Start >= xs2Start - clustdist) && (s2Start <= xs2End + clustdist) ||
                     (xs2Start >= s2Start - clustdist) && (xs2Start <= s2End + clustdist)) {
                    // refid needs to match
                    if((clusterlist[i].elements[0].refid == clusterlist[j].elements[0].refid) &&
                       (clusterlist[i].elements[0].refid == clusterlist[j].elements[0].refid)) {
                        // same with strand
                        if((clusterlist[i].elements[0].flag == clusterlist[j].elements[0].flag) &&
                           (clusterlist[i].elements[1].flag == clusterlist[j].elements[1].flag)) {
                            Cluster ncl = clusterlist[j];
                            ncl.count++;
                            // determine overlapping boundaries
                            if(s1Start < xs1Start) { ncl.elements[0].start = s1Start; }
                            if(s1End > xs1End) { ncl.elements[0].end = s1End; }
                            if(s2Start < xs2Start) { ncl.elements[1].start = s2Start; }
                            if(s2End > xs2End) { ncl.elements[1].end = s2End; }

                            clusterlist[j] = ncl;
                            clusterlist.erase(clusterlist.begin()+i); // remove cluster
                            --i;
                            break;
                        }
                    }
                }
            }
        }
    }
}


void Clustering::start(pt::ptree sample) {
    pt::ptree input = sample.get_child("input");
    std::string splits = input.get<std::string>("splits");

    std::cout << helper::getTime() << "Process sample: " << splits << "\n";
    iterate(splits);
}