#include "Analysis.hpp"

//
Analysis::Analysis(po::variables_map _params) : params(_params) {
    std::string line;
    std::ifstream anno;

    anno.open(params["features"].as<std::string>());
    // parse annotations
    if(!anno.is_open()) {
        perror("Error open");
        exit(EXIT_FAILURE);
    }
    while(getline(anno, line)) {
        if(line[0] == '#') { continue; }

        std::vector<std::string> tokens;
        std::istringstream iss(line);
        std::string token;
        while(std::getline(iss, token, '\t')) {
            tokens.push_back(token);
        }

        // boundaries
        std::pair<int,int> bnds = std::make_pair(std::stoi(tokens[3]),std::stoi(tokens[4]));
        std::pair<std::pair<int,int>,std::string> con = std::make_pair(bnds,"strand="+tokens[6]+";"+tokens[8]);

        std::size_t start_position;
        if(features.count(tokens[0]) == 0) { // check of RefName is in map
            if(start_position != std::string::npos) {
                features.insert(std::pair<std::string, std::vector<std::pair<std::pair<int,int>,std::string>>>(tokens[0],{con}));
            }
        } else {
            features[tokens[0]].push_back(con);
        }
    }

    // read file
    fs::path freqFile = fs::path(params["outdir"].as<std::string>()) / fs::path("frequency.txt");
    std::ifstream freqHandle;
    freqHandle.open(freqFile.string());

    // iterate and store into frquency map
    while(getline(freqHandle, line)) {
        std::vector<std::string> tokens;
        std::istringstream iss(line);
        std::string token;
        while(std::getline(iss, token, '\t')) {
            tokens.push_back(token);
        }
        frequency.insert(std::pair<std::string,int>(tokens[0],std::stoi(tokens[1])));
    }
}

//
void Analysis::createCountTable() {
    std::map<std::tuple<std::string,std::string,std::string,std::string>,std::vector<std::tuple<int,std::vector<float>,std::vector<float>>>> counts;

    std::ifstream intsfile;
    for(unsigned i=0;i<interPaths.size();++i) {
        intsfile.open(interPaths[i]);

        std::cout << interPaths[i] << std::endl;

        if(!intsfile.is_open()) {
            perror("Error open");
            exit(EXIT_FAILURE);
        }

        std::string line;
        while(getline(intsfile, line)){  //read data from file object and put it into string.
            if(line[0] == '#') { continue; }

            std::vector<std::string> tokens;
            std::istringstream iss(line);
            std::string token;
            while(std::getline(iss, token, '\t')) {
                tokens.push_back(token);
            }

            std::tuple<std::string,std::string,std::string,std::string> key1;
            std::tuple<std::string,std::string,std::string,std::string> key2;
            key1 = std::make_tuple(tokens[5],tokens[8],tokens[13],tokens[16]);
            key2 = std::make_tuple(tokens[5],tokens[8],tokens[13],tokens[16]);

            //
            if(counts.count(key1) == 0 && counts.count(key2) == 0) {
                std::vector<std::tuple<int,std::vector<float>,std::vector<float>>> content;

                for(unsigned j=0; j<interPaths.size();++j) {
                    std::tuple<int,std::vector<float>,std::vector<float>> vals{0,{},{}};
                    content.push_back(vals);
                }
                std::get<0>(content[i]) = 1;
                std::vector<float> en = std::get<1>(content[i]);
                std::vector<float> cm = std::get<2>(content[i]);

                en.push_back(std::stof(tokens[17]));
                cm.push_back(std::stof(tokens[18]));

                std::get<1>(content[i]) = en;
                std::get<2>(content[i]) = cm;

                counts.insert(std::pair<std::tuple<std::string,std::string,std::string,std::string>,std::vector<std::tuple<int,std::vector<float>,std::vector<float>>>>(key1,content));

            } else {
                if(counts.count(key1) > 0) {
                    std::tuple<int,std::vector<float>,std::vector<float>> oldVal;
                    oldVal = counts[key1][i];
                    std::get<0>(oldVal) = std::get<0>(oldVal) + 1;

                    std::vector<float> oldValEnergy = std::get<1>(oldVal);
                    std::vector<float> oldValCmpl = std::get<2>(oldVal);

                    oldValEnergy.push_back(std::stof(tokens[17]));
                    oldValCmpl.push_back(std::stof(tokens[18]));

                    std::get<1>(oldVal) = oldValEnergy;
                    std::get<2>(oldVal) = oldValCmpl;

                    counts[key1][i] = oldVal;
                } else {
                    if(counts.count(key2) > 0) {
                        std::tuple<int,std::vector<float>,std::vector<float>> oldVal;
                        oldVal = counts[key2][i];
                        std::get<0>(oldVal) = std::get<0>(oldVal) + 1;

                        std::vector<float> oldValEnergy = std::get<1>(oldVal);
                        std::vector<float> oldValCmpl = std::get<2>(oldVal);

                        oldValEnergy.push_back(std::stof(tokens[17]));
                        oldValCmpl.push_back(std::stof(tokens[18]));

                        std::get<1>(oldVal) = oldValEnergy;
                        std::get<2>(oldVal) = oldValCmpl;

                        counts[key2][i] = oldVal;
                    }
                }
            }
        }
    }


    // write back to file
    std::string outfile = params["outdir"].as<std::string>();
    fs::path outPath = fs::path(outfile) / "allints.txt";

    std::ofstream outFileHandle;
    outFileHandle.open(outPath.string());

    outFileHandle << "RNA1\tRNA2\tRNA1orientation\tRNA2orientation\t";

    for(unsigned i=0;i<interPaths.size();i++) {
        outFileHandle << fs::path(interPaths[i]).stem().string() << "_counts\t";
        outFileHandle << fs::path(interPaths[i]).stem().string() << "_ges\t";
        outFileHandle << fs::path(interPaths[i]).stem().string() << "_gcs\t";
    }
    outFileHandle << "\n";

    for(auto it=counts.begin(); it!=counts.end(); ++it) {
        outFileHandle << std::get<0>(it->first) << "\t";
        outFileHandle << std::get<2>(it->first) << "\t";
        outFileHandle << std::get<1>(it->first) << "\t";
        outFileHandle << std::get<3>(it->first) << "\t";

        //
        for(unsigned z=0;z<it->second.size();++z) {
            outFileHandle << std::get<0>(it->second[z]) << "\t";

            std::vector<float> vecNrg = std::get<1>(it->second[z]);
            std::vector<float> vecCpl = std::get<2>(it->second[z]);

            // determine ges
            float ges = 0.0;
            if(vecNrg.size() == 1) {
                ges = vecNrg[0];
            } else {
                if(vecNrg.size() > 1) {
                    auto m = vecNrg.begin() + vecNrg.size()/2;
                    std::nth_element(vecNrg.begin(), m, vecNrg.end());
                    double min = *min_element(vecNrg.begin(), vecNrg.end());
                    ges = vecNrg[vecNrg.size()/2];
                }
            }
            outFileHandle << ges << "\t";

            // determine gcs
            float gcs = 0.0;
            if(vecCpl.size() == 1) {
                gcs = vecCpl[0];
            } else {
                if(vecCpl.size() > 1) {
                    auto m = vecCpl.begin() + vecCpl.size()/2;
                    std::nth_element(vecCpl.begin(), m, vecCpl.end());
                    double min = *min_element(vecCpl.begin(), vecCpl.end());
                    gcs = vecCpl[vecCpl.size()/2];
                }
            }
            outFileHandle << gcs << "\t";
        }
        outFileHandle << "\n";
    }
    outFileHandle.close();
}



//
std::string Analysis::retrieveTagValue(std::string tags, std::string tagName, std::string oldValue) {
    std::size_t start_position = tags.find(tagName+"=");
    // gene name
    if(start_position != std::string::npos) {
        std::string sub = tags.substr(start_position+tagName.size()+1,tags.length());
        std::size_t end_position = sub.find(";");
        oldValue = sub.substr(0,end_position);
    }
    return oldValue;
}


float Analysis::calc_pvalue(int x, int n, float p) {
    // simulate draw from binomial distribution
    std::binomial_distribution<int> binomial(n, p);
    // assign p value
    float pval = 1.0 - cdf(binomial, x);
    return pval;
}

void Analysis::start(pt::ptree sample) {
    // retrieve input and output files
    pt::ptree input = sample.get_child("input");
    std::string splits = input.get<std::string>("splits");
    pt::ptree output = sample.get_child("output");
    std::string interactions = output.get<std::string>("interactions");


    interPaths.push_back(interactions);

    // input .sam record
    seqan3::alignment_file_input fin{splits,
                                     seqan3::fields<seqan3::field::id,
                                             seqan3::field::flag,
                                             seqan3::field::ref_id,
                                             seqan3::field::ref_offset,
                                             seqan3::field::seq,
                                             seqan3::field::tags>{}};

    std::vector<seqan3::sam_flag> flags;
    std::vector<std::string> refIDs;
    std::vector<std::optional<int32_t>> refOffsets;

    std::vector<size_t> ref_lengths{};
    for(auto &info : fin.header().ref_id_info) {
        ref_lengths.push_back(std::get<0>(info));
    }
    std::deque<std::string> ref_ids = fin.header().ref_ids();

//    uint32_t flag, start, end;

    // open file
    std::ofstream outInts;
    outInts.open(interactions);

    std::string entry; // stores output to write to file

    // variables for interactions file
    std::string qNAME, flag;
    uint32_t start, end;
    std::string geneID, geneName, product, annoStrand;
    float hybnrg, cmpl;

    seqan3::sam_tag_dictionary tags;

    outInts << "#QNAME\tSegment1Strand\tSegment1Start\tSegment1End\tSegment1RefName"
    outInts << "\tSegment1Name\tSegment1AnnoStrand\tSegment1Product\tSegment1Orientation"
    outInts << "\tSegment2Strand\tSegment2Start\tSegment2End\tSegment2RefName\tSegment2Name"
    outInts << "\tSegment2AnnoStrand\tSegment2Product\tSegment1Orientation\tenergy\tcomplementarity\n";

    int segCnt = 0;
    int segCntMatch = 0;
    int segFound = 0;

    for(auto && rec : fin | seqan3::views::chunk(2)) {
        entry = "";
        hybnrg = 0.0;
        cmpl = 0.0;
        for(auto & split : rec) {
            //           seqan3::debug_stream << split << std::endl;

            qNAME = seqan3::get<seqan3::field::id>(split);
            flag = "+";
            if(static_cast<bool>(seqan3::get<seqan3::field::flag>(split)
                                 & seqan3::sam_flag::on_reverse_strand)) {
                flag = "-";
            }

            // start & end
            start = seqan3::get<seqan3::field::ref_offset>(split).value();
            end = start + seqan3::get<seqan3::field::seq>(split).size()-1;

            // refID
            std::optional<int32_t> refIDidx = seqan3::get<seqan3::field::ref_id>(split);
            std::string refID = ref_ids[refIDidx.value()];

            tags = seqan3::get<seqan3::field::tags>(split);
            auto nrg = tags["XE"_tag];
            auto cpl = tags["XC"_tag];
            hybnrg = std::get<float>(nrg);
            cmpl = std::get<float>(cpl);

            std::vector<std::pair<std::pair<int,int>,std::string>> feat;
            if(features.count(refID) > 0) {
                feat = features[refID];

                std::pair<std::string, std::string> geneNames;
                geneID = ".";
                geneName = ".";
                product = ".";
                annoStrand = ".";
                for(unsigned i=0;i<feat.size();++i) {
                    if((start >= feat[i].first.first && start <= feat[i].first.second) ||
                       (end >= feat[i].first.first && end <= feat[i].first.second)) {

                        // make sure that overlap is in exon
                        if(params["splicing"].as<std::bitset<1>() && feat[i].second.find("exon") == std::string::npos) {
                            continue;
                        }
                        geneID = retrieveTagValue(feat[i].second, "ID", geneName);
                        geneName = retrieveTagValue(feat[i].second, "gene", geneName);
                        product = retrieveTagValue(feat[i].second, "product", product);
                        annoStrand = retrieveTagValue(feat[i].second, "strand", annoStrand);

                        segFound = 1; // found a match (in annotation) for segment
                    }
                }

                // search for geneNames in frequency map
                geneNames = std::make_pair(retrieveTagValue(feat[0].second, "gene", geneName), retrieveTagValue(feat[1].second, "gene", geneName));

                float value;
                // check if geneNames first in frequency map
                if(frequency.count(geneNames.first) > 0 && frequency.count(geneNames.second) > 0) {
                    if (geneNames.first == geneNames.second) {
                        value = frequency[geneNames.first] * frequency[geneNames.second];
                    } else {
                        value = 2 * (frequency[geneNames.first] * frequency[geneNames.second]);
                    }
                } else{
                    value = 0;
                }

                // calculate p-value
                int draws = frequency[geneNames.first] + frequency[geneNames.second];
                float pval = calc_pvalue(segCntMatch, draws, value);


                if(segFound == 1) {
                    if(segCnt == 0) { entry += qNAME + "\t";}
                    entry += flag + "\t";
                    entry += std::to_string(start) + "\t";
                    entry += std::to_string(end) + "\t";
                    entry += refID + "\t";
                    entry += geneName + "\t";
                    entry += annoStrand + "\t";
                    entry += "\"" + product  + "\"\t";

                    if(flag == annoStrand) {
                        entry += "sense\t";
                    } else {
                        entry += "antisense\t";
                    }

                    entry += pval + "\t";

                    segCntMatch++;
                    segFound = 0;
                }
            }
            segCnt++; // on to the next segment
        }

        if(segCntMatch == 2) {
            outInts << entry << hybnrg << "\t" << cmpl << "\n";
        }

        segCnt = 0; // reset segment
        segCntMatch = 0;
    }

    outInts.close();
}


