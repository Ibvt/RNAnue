#include <iostream>
#include "Clustering.hpp"

Clustering::Clustering(po::variables_map params) : 
    params(params) {
        result = {}; 
}


void Clustering::iterate(std::string splits, std::string clusters) {

    // input .sam record
    seqan3::alignment_file_input fin{splits, 
		seqan3::fields<seqan3::field::id,
		                seqan3::field::flag,
                        seqan3::field::ref_id,
						seqan3::field::ref_offset,
						seqan3::field::seq,
						seqan3::field::tags>{}};

    int chunks = 5;
    std::vector<Cluster> subset;

    //
    std::vector<seqan3::sam_flag> flags;
    std::vector<std::string> refIDs;
    std::vector<std::optional<int32_t>> refOffsets;
    
    std::vector<size_t> ref_lengths{};
	for(auto &info : fin.header().ref_id_info) {
		ref_lengths.push_back(std::get<0>(info));
	}
    std::deque<std::string> ref_ids = fin.header().ref_ids();

    uint32_t flag, start, end;
    for(auto && rec : fin | seqan3::views::chunk(2)) {
        Cluster cl;
        for(auto & split : rec) {
            std::optional<int32_t> refID = seqan3::get<seqan3::field::ref_id>(split); 
            uint32_t flag{0}; // SAMFLAG
            if(static_cast<bool>(seqan3::get<seqan3::field::flag>(split) & seqan3::sam_flag::on_reverse_strand)) {
                flag = 16;
            } 
            
            // start & end
            uint32_t start = seqan3::get<seqan3::field::ref_offset>(split).value();
            uint32_t end = start + seqan3::get<seqan3::field::seq>(split).size()-1;

            Segment s(ref_ids[refID.value()], flag, start, end);
            cl.elements.push_back(s);
        }

        // always have the left segment occurring first
        if(cl.elements[0].start > cl.elements[1].start) {
            Segment cl0 = cl.elements[0];
            cl.elements[0] = cl.elements[1];
            cl.elements[1] = cl0;
        }
        subset.push_back(cl);

        std::cout << "size of subset: " << subset.size() << std::endl;

        if(subset.size() == chunks) {
            overlaps(subset);
            result.insert(result.end(),subset.begin(),subset.end());
            std::cout << "results size: " << result.size() << std::endl;
            subset.clear();
        }
    }
    if(subset.size() > 0) {
        overlaps(subset);
        result.insert(result.end(),subset.begin(),subset.end());
        std::cout << "results size: " << result.size() << std::endl;
        subset.clear();
    }

    overlaps(result);
    std::cout << "final results size: " << result.size() << std::endl;
}


//
bool Clustering::startPosCmp(Cluster &a, Cluster &b) {
    return a.elements[0].start < b.elements[0].start;
}


void Clustering::sumup() {
    std::cout << "sumup" << std::endl;
    std::cout << result.size() << std::endl;
    /*
    for(unsigned i=0;clusters.size();++i) {
        //std::cout << clusters[i].count << std::endl;

    }*/

}
void Clustering::overlaps(std::vector<Cluster> &clusterlist) {
    std::sort(clusterlist.begin(), clusterlist.end());

    uint32_t s1Start, s1End, s2Start, s2End;
    uint32_t xs1Start, xs1End, xs2Start, xs2End;

    for(unsigned i=0;i<clusterlist.size();++i) {
//        std::cout << clusterlist[i].elements[0].start << std::endl;
        for(unsigned j=i+1; j<clusterlist.size();++j) {
            
            uint32_t s1Start = clusterlist[i].elements[0].start+1;
            uint32_t s1End = clusterlist[i].elements[0].end+1;
            uint32_t s2Start = clusterlist[i].elements[1].start+1;
            uint32_t s2End = clusterlist[i].elements[1].end+1;
            
            uint32_t xs1Start = clusterlist[j].elements[0].start+1;
            uint32_t xs1End = clusterlist[j].elements[0].end+1;
            uint32_t xs2Start = clusterlist[j].elements[1].start+1;
            uint32_t xs2End = clusterlist[j].elements[1].end+1;

            /*           
            std::cout << "######### i " << std::endl;        
            std::cout << "first" << std::endl;
            std::cout << "start: " << s1Start << std::endl;
            std::cout << "end: " << s1End << std::endl;
            std::cout << "second" << std::endl;
            std::cout << "start: " << s2Start << std::endl;
            std::cout << "end: " << s2End << std::endl;

            std::cout << "######### j " << std::endl;        
            std::cout << "first" << std::endl;
            std::cout << "start: " << xs1Start << std::endl;
            std::cout << "end: " << xs1End << std::endl;
            std::cout << "second" << std::endl;
            std::cout << "start: " << xs2Start << std::endl;
            std::cout << "end: " << xs2End << std::endl;
            */

            if((xs1Start <= s1End)) { // first split matches
                if((s2Start >= xs2Start) && (s2Start <= xs2End) ||
                        (xs2Start >= s2Start) && (xs2Start <= s2End)) {

                    // refid needs to match
                    if((clusterlist[i].elements[0].refid == clusterlist[j].elements[0].refid) &&
                            (clusterlist[i].elements[0].refid == clusterlist[j].elements[0].refid)) {

                        // .. same with strand
                        if((clusterlist[i].elements[0].flag == clusterlist[j].elements[0].flag) &&
                                (clusterlist[i].elements[1].flag == clusterlist[j].elements[1].flag)) {

                            Cluster ncl = clusterlist[j];
                            ncl.count += 1;

                            if(s1Start < xs1Start) {
                                ncl.elements[0].start = s1Start;
                            }
                            if(s1End > xs1End) {
                                ncl.elements[0].end = s1End;
                            }
                            if(s2Start < xs2Start) {
                                ncl.elements[1].start = s2Start;
                            }
                            if(s2End > xs2End) {
                                ncl.elements[1].end = s2End;
                            }

                            clusterlist[j] = ncl;

                            // remove cluster
                            clusterlist.erase(clusterlist.begin()+i);
                            //std::remove(clusterlist.begin(),clusterlist.end(),clusterlist[i]);


                        }
                    }
                }
            }
        }
    }
}



void Clustering::start(pt::ptree sample) {
    pt::ptree input = sample.get_child("input");
    pt::ptree output = sample.get_child("output");

    std::string splits = input.get<std::string>("splits");
    std::string clusters = output.get<std::string>("clusters");

    std::cout << splits << std::endl;

    iterate(splits,clusters);

    //std::cout << "within start" << splits << std::endl;
}
