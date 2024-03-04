#include "SplitReadCalling.hpp"


SplitReadCalling::SplitReadCalling(po::variables_map params) :
    params(params) {
    readscount = 0;
    alignedcount = 0;
    splitscount = 0;
    msplitscount = 0;
    nsurvivedcount = 0;

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
}


~SplitReadCalling::SplitReadCalling() {
    // output file
    std::ofstream ffile;
    fs::path freqfile = fs::path(params["outdir"].as<std::string>()) / fs::path("frequency.txt");

    // iterate through map
    for(auto &it : frequency) {
        if(ffile.is_open()) {
            float freq = (float)it.second / (float)readscount;
            ffile << freq << "\t" << it.first << "\t" << it.second << std::endl;
        }
    }
    ffile.close();
}


//
void SplitReadCalling::iterate(std::string matched, std::string splits, std::string multsplits) {
    // input .sam record
    seqan3::alignment_file_input fin{matched,
                                     seqan3::fields<seqan3::field::id, seqan3::field::flag, seqan3::field::ref_id,
                                             seqan3::field::ref_offset, seqan3::field::mapq, seqan3::field::cigar,
                                             seqan3::field::seq, seqan3::field::tags, seqan3::field::alignment>{}};


    std::vector<size_t> ref_lengths{};
    for(auto &info : fin.header().ref_id_info) {
        ref_lengths.push_back(std::get<0>(info));
    }
    std::deque<std::string> ref_ids = fin.header().ref_ids();

    // output file for splits
    seqan3::alignment_file_output splitsfile{splits, ref_ids, ref_lengths,
                                             seqan3::fields<seqan3::field::id, seqan3::field::flag, seqan3::field::ref_id,
                                                     seqan3::field::ref_offset, seqan3::field::mapq, seqan3::field::cigar,
                                                     seqan3::field::seq, seqan3::field::tags>{}};


    // output file for multi splits
    seqan3::alignment_file_output multsplitsfile{multsplits, ref_ids, ref_lengths,
                                                 seqan3::fields<
                                                         seqan3::field::id,
                                                         seqan3::field::flag,
                                                         seqan3::field::ref_id,
                                                         seqan3::field::ref_offset,
                                                         seqan3::field::mapq,
                                                         seqan3::field::cigar,
                                                         seqan3::field::seq,
                                                         seqan3::field::tags>{}};


    std::string currentQNAME = "";
    //std::vector<std::string> splitrecord; // all split combinations
    using record_type = typename decltype(fin)::record_type;
    std::list<record_type> readrecords{};
    std::vector<std::list<record_type>> tasks;

    for (auto & rec : fin) {
        // ignore reads if unmapped and do not suffice length requirements
        if(static_cast<bool>(seqan3::get<seqan3::field::flag>(rec) & seqan3::sam_flag::unmapped)) {
            continue;
        }
        if(seqan3::get<seqan3::field::seq>(rec).size() <= params["minlen"].as<int>()) {
            continue;
        }

        std::string QNAME = seqan3::get<seqan3::field::id>(rec);
        if((currentQNAME != "") && (currentQNAME != QNAME)) {
            readscount++;
            process(readrecords, splitsfile, multsplitsfile);
            readrecords.clear();
            std::list<record_type>().swap(readrecords);

            // add for next iteration
            readrecords.push_back(rec);
            currentQNAME = QNAME;
            if(readscount % 100000 == 0) {
                progress(std::cout);
            }
        } else {
            readrecords.push_back(rec);
            currentQNAME = QNAME;
        }
    }
    if(!readrecords.empty()) {
        readscount++;
        process(readrecords, splitsfile, multsplitsfile);
        readrecords.clear();
        std::list<record_type>().swap(readrecords);
    }
    progress(std::cout);
}

// distribute calculations across cores
void SplitReadCalling::distribute(auto &tasks, auto &splits, auto &multsplits) {
    // distribute tasks
    unsigned i;
# pragma omp parallel for shared(splits, multsplits) num_threads(params["threads"].as<int>())
    for(i=0;i<tasks.size();++i) {
        process(tasks[i],splits,multsplits);
    }
}


void SplitReadCalling::process(auto &splitrecords, auto &splitsfile, auto &multsplitsfile) {
    Splts splits;

    seqan3::cigar_op cigarOp; // operation of cigar string
    uint32_t cigarOpSize;

    uint32_t startPosRead; uint32_t endPosRead; // absolute position in read/alignment (e.g., 1 to end)
    uint32_t startPosSplit; uint32_t endPosSplit; // position in split read (e.g., XX:i to XY:i / 14 to 20)

    seqan3::sam_tag_dictionary tags;

    int segmentNr = 0; // intialize active segment
    int splitID = 0; //

    std::string qname = "";
    std::vector<SamRecord> curated;

    // create object of putative splits
    std::map<int, std::vector<SamRecord>> putative;
    std::vector<std::pair<SamRecord,SamRecord>> splitSegments;

    //
    auto it = splitrecords.begin();
    while(it != splitrecords.end()) {
        // extract information
        qname = seqan3::get<seqan3::field::id>(*it); // QNAME
        seqan3::sam_flag flag = seqan3::get<seqan3::field::flag>(*it); // SAMFLAG
        std::optional<int32_t> refID = seqan3::get<seqan3::field::ref_id>(*it);
        std::optional<int32_t> refOffset = seqan3::get<seqan3::field::ref_offset>(*it);
        std::optional<uint8_t> qual = seqan3::get<seqan3::field::mapq>(*it);

        // CIGAR string
        std::vector<seqan3::cigar> cigar{seqan3::get<seqan3::field::cigar>(*it)};
        std::vector<seqan3::cigar> cigarSplit{}; // individual cigar for each split
        std::span<seqan3::dna5> seq = seqan3::get<seqan3::field::seq>(*it);

        // extract the tags (information about split reads)
        tags = seqan3::get<seqan3::field::tags>(*it);
        auto xhtag = tags.get<"XH"_tag>();
        auto xjtag = tags.get<"XJ"_tag>();
        auto xxtag = tags.get<"XX"_tag>();
        auto xytag = tags.get<"XY"_tag>();
//        auto xctag = tags.get<"XC"_tag>();

        if(xjtag >= 2) { // splits in SAMfile need to consists of at least 2 segments
            startPosRead = 1; // intialise start & end of reads
            endPosRead = 0;

            startPosSplit = xxtag; // determine the start position within split
            endPosSplit = xxtag-1;

            // check the cigar string for splits
            for(auto &cig : cigar) {
                // determine size and operator of cigar element
                cigarOpSize = get<uint32_t>(cig);
                cigarOp = get<seqan3::cigar_op>(cig);

                //
                if(cigarOp == 'N'_cigar_op) {
                    auto subSeq = seq | seqan3::views::slice(startPosRead-1,endPosRead);

                    // add properties to tags - lightweight
                    seqan3::sam_tag_dictionary ntags;
                    ntags.get<"XX"_tag>() = startPosSplit;
                    ntags.get<"XY"_tag>() = endPosSplit;
                    ntags["XN"_tag] = splitID;

                    filterSegments(*it, refOffset, cigarSplit, subSeq, ntags, curated);

                    // settings for prepare for next split
                    startPosSplit = endPosSplit+1;
                    startPosRead = endPosRead+1;
                    cigarSplit = {}; // new split - new CIGAR
                    refOffset.value() += cigarOpSize + endPosRead+1; // adjust leftmost mapping position

                    // change for next iteration
                    segmentNr++; // counts as segment of split

                } else {
                    // deletion does not account for length in split
                    seqan3::cigar cigarElement{cigarOpSize, cigarOp};
                    //seqan3::debug_stream << "cigarElement: " << cigarElement << std::endl;

                    // exclude soft clipping from alignment
                    if(cigarOp == 'S'_cigar_op && params["exclclipping"].as<std::bitset<1>>() == 1) {
                        if(cigarSplit.size() == 0) {
                            startPosRead += cigarOpSize;
                            endPosRead += cigarOpSize;
                            startPosSplit += cigarOpSize;
                            endPosSplit += cigarOpSize;
                        }
                    } else {
                        if(cigarOp != 'D'_cigar_op) {
                            endPosRead += cigarOpSize;
                            endPosSplit += cigarOpSize;
                        }
                        cigarSplit.push_back(cigarElement);
                    }
                }
            }

            std::vector<std::pair<std::pair<int,int>,std::string>> feat;
            if(features.count(refID) > 0) {
                feat = features[refID];

                geneName = ".";
                product = ".";
                annoStrand = ".";
                for(unsigned i=0;i<feat.size();++i) {
                    if((start >= feat[i].first.first && start <= feat[i].first.second) ||
                       (end >= feat[i].first.first && end <= feat[i].first.second)) {

                        geneName = retrieveTagValue(feat[i].second, "ID", geneName);
                        geneName = retrieveTagValue(feat[i].second, "gene", geneName);
                        product = retrieveTagValue(feat[i].second, "product", product);
                        annoStrand = retrieveTagValue(feat[i].second, "strand", annoStrand);

                        segFound = 1; // found a match (in annotation) for segment

                        //
                        if (it != frequency.end()) {
                            // Key exists, increment value by 1
                            it->second++;
                        }
                    }
                }

            auto subSeq = seq | seqan3::views::slice(startPosRead-1,endPosRead);

            seqan3::sam_tag_dictionary ntags;
            ntags.get<"XX"_tag>() = startPosSplit;
            ntags.get<"XY"_tag>() = endPosSplit;
            ntags["XN"_tag] = splitID;

            filterSegments(*it, refOffset, cigarSplit, subSeq, ntags, curated);
            segmentNr++;

            if(segmentNr == xjtag) {
                std::map<int, std::vector<SamRecord>>::iterator itPutative = putative.begin();
                putative.insert(itPutative, std::make_pair(splitID,curated));
                curated.clear();
                segmentNr = 0;
                ++splitID;
            }
            ++it;
        } else {
            ++it;
        }
    }

    if(putative.size() > 0) { // split reads found
        std::vector<SamRecord> p1;
        std::vector<SamRecord> p2;

        seqan3::sam_tag_dictionary p1Tags;
        seqan3::sam_tag_dictionary p2Tags;

        // value for complementarity and hybridization energy
        std::pair<double,double> filters;

        std::string p1Qname;
        std::string p2Qname;

        std::map<int, std::vector<SamRecord>>::iterator itSplits = putative.begin();
        for(itSplits;itSplits != putative.end(); ++itSplits) {
            std::vector<SamRecord> splits = itSplits->second;
            if(splits.size() > 1) { // splits consists at least of 2 segments
                for(unsigned i=0;i<splits.size();++i) {
                    for(unsigned j=i+1;j<splits.size();++j) {
                        p1Tags = seqan3::get<seqan3::field::tags>(splits[i]);
                        p2Tags = seqan3::get<seqan3::field::tags>(splits[j]);

                        auto p1Start = p1Tags.get<"XX"_tag>();
                        auto p1End = p1Tags.get<"XY"_tag>();
                        auto p2Start = p2Tags.get<"XX"_tag>();
                        auto p2End = p2Tags.get<"XY"_tag>();

                        // prevent an overlap between reads position -> same segment of read
                        if((p2Start >= p1Start && p2Start <= p1End) || (p1Start >= p2Start && p1Start <= p2End)) {
                            continue;
                        }

                        if(p1.empty() && p2.empty()) {
                            TracebackResult initCmpl = complementarity(seqan3::get<seqan3::field::seq>(splits[i]), seqan3::get<seqan3::field::seq>(splits[j]));
                            double initHyb = hybridize(seqan3::get<seqan3::field::seq>(splits[i]), seqan3::get<seqan3::field::seq>(splits[j]));
                            if(!initCmpl.a.empty()) { // split read survived complementarity & sitelenratio cutoff
                                if(initHyb <= params["nrgmax"].as<double>()) {
                                    filters = std::make_pair(initCmpl.cmpl, initHyb);
                                    addComplementarityToSamRecord(splits[i], splits[j], initCmpl);
                                    addHybEnergyToSamRecord(splits[i], splits[j], initHyb);
                                    splitSegments.push_back(std::make_pair(splits[i],splits[j]));
                                }
                            }
                        } else {
                            TracebackResult cmpl = complementarity(seqan3::get<seqan3::field::seq>(splits[i]), seqan3::get<seqan3::field::seq>(splits[j]));
                            double hyb = hybridize(seqan3::get<seqan3::field::seq>(splits[i]), seqan3::get<seqan3::field::seq>(splits[j]));

                            if(!cmpl.a.empty()) { // check that there is data for complementarity / passed cutoffs
                                if(cmpl.cmpl > filters.first) { // complementarity is higher
                                    splitSegments.clear();
                                    filters = std::make_pair(cmpl.cmpl, hyb);
                                    addComplementarityToSamRecord(splits[i], splits[j], cmpl);
                                    addHybEnergyToSamRecord(splits[i], splits[j], hyb);
                                    splitSegments.push_back(std::make_pair(splits[i],splits[j]));
                                } else{
                                    if(cmpl.cmpl == filters.first) {
                                        if(hyb > filters.second) {
                                            splitSegments.clear();
                                            filters = std::make_pair(cmpl.cmpl, hyb);
                                            addComplementarityToSamRecord(splits[i], splits[j], cmpl);
                                            addHybEnergyToSamRecord(splits[i], splits[j], hyb);
                                            splitSegments.push_back(std::make_pair(splits[i],splits[j]));
                                        }
                                    } else {
                                        if(hyb == filters.second) { // same cmpl & hyb -> ambigious read!
                                            addComplementarityToSamRecord(splits[i], splits[j], cmpl);
                                            addHybEnergyToSamRecord(splits[i], splits[j], hyb);
                                            splitSegments.push_back(std::make_pair(splits[i],splits[j]));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    {
        // write to file
        if(splitSegments.size() == 1) {
            writeSamFile(splitsfile, splitSegments);
            splitscount++;
        } else {
            if(splitSegments.size() > 1) {
                writeSamFile(multsplitsfile, splitSegments);
                msplitscount++;
            }
        }
        std::vector<std::pair<SamRecord,SamRecord>>().swap(splitSegments);
    }
}


// write splits back to file
void SplitReadCalling::writeSamFile(auto &samfile, std::vector<std::pair<SamRecord,SamRecord>> &splits) {
    for(unsigned i=0;i<splits.size();++i) {
        samfile.emplace_back(
                seqan3::get<seqan3::field::id>(splits[i].first),
                seqan3::get<seqan3::field::flag>(splits[i].first),
                seqan3::get<seqan3::field::ref_id>(splits[i].first),
                seqan3::get<seqan3::field::ref_offset>(splits[i].first),
                seqan3::get<seqan3::field::mapq>(splits[i].first).value(),
                seqan3::get<seqan3::field::cigar>(splits[i].first),
                seqan3::get<seqan3::field::seq>(splits[i].first),
                seqan3::get<seqan3::field::tags>(splits[i].first));

        samfile.emplace_back(
                seqan3::get<seqan3::field::id>(splits[i].second),
                seqan3::get<seqan3::field::flag>(splits[i].second),
                seqan3::get<seqan3::field::ref_id>(splits[i].second),
                seqan3::get<seqan3::field::ref_offset>(splits[i].second),
                seqan3::get<seqan3::field::mapq>(splits[i].second).value(),
                seqan3::get<seqan3::field::cigar>(splits[i].second),
                seqan3::get<seqan3::field::seq>(splits[i].second),
                seqan3::get<seqan3::field::tags>(splits[i].second));
    }
}


void SplitReadCalling::filterSegments(auto &splitrecord, std::optional<int32_t> &refOffset,
                                      std::vector<seqan3::cigar> &cigar, std::span<seqan3::dna5> &seq,
                                      seqan3::sam_tag_dictionary &tags, std::vector<SamRecord> &curated) {
    SamRecord segment{}; // create new SamRecord

    //
    seqan3::get<seqan3::field::id>(segment) = seqan3::get<seqan3::field::id>(splitrecord);
    seqan3::sam_flag flag{0}; // SAMFLAG
    if(static_cast<bool>(seqan3::get<seqan3::field::flag>(splitrecord) & seqan3::sam_flag::on_reverse_strand)) {
        flag = seqan3::sam_flag{16};
    }

    seqan3::get<seqan3::field::flag>(segment) = flag;
    seqan3::get<seqan3::field::ref_id>(segment) = seqan3::get<seqan3::field::ref_id>(splitrecord);
    seqan3::get<seqan3::field::ref_offset>(segment) = refOffset;
    seqan3::get<seqan3::field::mapq>(segment) = seqan3::get<seqan3::field::mapq>(splitrecord);
    seqan3::get<seqan3::field::cigar>(segment) = cigar;
    seqan3::get<seqan3::field::seq>(segment) = seq;
    seqan3::get<seqan3::field::tags>(segment) = tags;

    if(seq.size() <= params["minfraglen"].as<int>()) {
        return;
    }

    // std::cout << "length " << seq.size() << std::endl;
    curated.push_back(segment);
}


void SplitReadCalling::addFilterToSamRecord(SamRecord &rec, std::pair<float,float> filters) {
    seqan3::get<seqan3::field::tags>(rec)["XC"_tag] = filters.first;
    seqan3::get<seqan3::field::tags>(rec)["XE"_tag] = filters.second;
}

void SplitReadCalling::addComplementarityToSamRecord(SamRecord &rec1, SamRecord &rec2, TracebackResult &res) {
    // number of matches
    seqan3::get<seqan3::field::tags>(rec1)["XM"_tag] = res.matches;
    seqan3::get<seqan3::field::tags>(rec2)["XM"_tag] = res.matches;

    // length of alignment
    seqan3::get<seqan3::field::tags>(rec1)["XL"_tag] = res.length;
    seqan3::get<seqan3::field::tags>(rec2)["XL"_tag] = res.length;

    // complementarity
    seqan3::get<seqan3::field::tags>(rec1)["XC"_tag] = (float)res.cmpl;
    seqan3::get<seqan3::field::tags>(rec2)["XC"_tag] = (float)res.cmpl;

    // sitelenratio
    seqan3::get<seqan3::field::tags>(rec1)["XR"_tag] = (float)res.ratio;
    seqan3::get<seqan3::field::tags>(rec2)["XR"_tag] = (float)res.ratio;

    // alignments
    seqan3::get<seqan3::field::tags>(rec1)["XA"_tag] = res.a;
    seqan3::get<seqan3::field::tags>(rec2)["XA"_tag] = res.b;

    // score
    seqan3::get<seqan3::field::tags>(rec1)["XS"_tag] = res.score;
    seqan3::get<seqan3::field::tags>(rec2)["XS"_tag] = res.score;
}

//
void SplitReadCalling::addHybEnergyToSamRecord(SamRecord &rec1, SamRecord &rec2, double &hyb) {
    seqan3::get<seqan3::field::tags>(rec1)["XE"_tag] = (float)hyb;
    seqan3::get<seqan3::field::tags>(rec2)["XE"_tag] = (float)hyb;
}

//
TracebackResult SplitReadCalling::complementarity(std::span<seqan3::dna5> &seq1, std::span<seqan3::dna5> &seq2) {
    std::string seq1_str;
    std::string seq2_str;

    if(seq1.size() == 0 || seq2.size() == 0) {
        return {"","",0,0,0,-1.0,0.0};
    }

    // rna1 (reverse)
    for(unsigned z=seq1.size();z-- > 0;) {
        seq1_str += seq1[z].to_char();
    }

    // rna2
    for(unsigned y=0;y<seq2.size();++y) {
        seq2_str += seq2[y].to_char();
    }

    ScoringMatrix matrix = createScoringMatrix(seq1_str.c_str(),seq2_str.c_str(), 1, -1, 2);
    std::vector<TracebackResult> res = traceback(matrix, seq1_str.c_str(),seq2_str.c_str());

    int idx = -1;
    int cmpl = -1;
    // filter (multiple) alignments
    if(res.size() > 1) {
        for(unsigned i=0;i<res.size();++i) {
            // check if sitelenratio & complementarity exceed cutoffs
            if(params["sitelenratio"].as<double>() <= res[i].ratio &&
               params["cmplmin"].as<double>() <= res[i].cmpl) {
                if(res[i].cmpl > cmpl) {
                    idx = i;
                }
            }
        }
    }

//    printMatrix(matrix); // print matrix
    freeMatrix(&matrix);
    if(idx != -1) {
        return res[idx];
    } else {
        return {"","",0,0,0,0.0,0.0};
    }
}


double SplitReadCalling::hybridize(std::span<seqan3::dna5> &seq1, std::span<seqan3::dna5> &seq2) {
    std::string rna1str = "";
    std::string rna2str = "";

    for(unsigned i=0;i<seq1.size();++i) { rna1str += seq1[i].to_char(); }
    for(unsigned i=0;i<seq2.size();++i) { rna2str += seq2[i].to_char(); }

    std::string hyb = "echo '" + rna1str + "&" + rna2str + "' | RNAcofold";

    const char* call = hyb.c_str();
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(call, "r"), pclose);

    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }

    std::stringstream ss(result);
    std::vector<std::string> tokens;
    std::string tmp;

    while(getline(ss,tmp,'\n')) {
        tokens.push_back(tmp);
    }

    std::string dotbracket = tokens[1].substr(0,tokens[1].find(' '));
    std::regex regexp("-?[[:digit:]]{1,2}\\.[[:digit:]]{1,2}");
    std::smatch matches;

    std::regex_search(tokens[1], matches, regexp);
    return stod(matches[0]);
}


// calculate rows in SAMfile (exclusing header)
int SplitReadCalling::countSamEntries(std::string file, std::string command) {
    std::string entries = "awk '!/^@/ {print $0}' " + file + " " + command + " | wc -l";
    const char* call = entries.c_str();
    std::array<char, 128> buffer;

    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(call, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return std::stoi(result);
}



// adds suffix to filename (optional:
std::string SplitReadCalling::addSuffix(std::string _file, std::string _suffix, std::vector<std::string> _keywords) {
    int keyPos, tmpKeyPos = -1; // buffer the positions of the keywords
    int dotPos = _file.find("."); // determine position of the dot
    keyPos = dotPos;
    if(!_keywords.empty()) {
        for(unsigned i=0;i<_keywords.size();++i) {
            tmpKeyPos = _file.find(_keywords[i]);
            if(tmpKeyPos != -1) { // key could be found
                keyPos = tmpKeyPos;
            }
        }
    }
    std::string newFile = _file.substr(0,keyPos);
    newFile += _suffix;
    newFile += _file.substr(dotPos,_file.size());

    return newFile;
}

void SplitReadCalling::createDir(fs::path path) {
    if(fs::exists(path)) {
        fs::remove_all(path);
    }
    fs::create_directory(path);
}

void SplitReadCalling::progress(std::ostream& out) {
    auto ctime = std::chrono::system_clock::now(); // get current time
    std::time_t end_time = std::chrono::system_clock::to_time_t(ctime);
    out << "\t\t" << std::ctime(&end_time);
    out << "\t\t" << readscount << " reads processed\n";
}


//
std::vector<std::vector<fs::path>> SplitReadCalling::splitInputFile(std::string matched, std::string splits, int entries) {
    fs::path matchedFilePath = fs::path(matched);
    fs::path splitsFilePath = fs::path(splits);

    // create tmp folder (for inputs)
    fs::path tmpIn = matchedFilePath.parent_path() / "tmpIn";
    fs::path tmpSplit = splitsFilePath.parent_path() / "tmpSplit";
    fs::path tmpMsplit = splitsFilePath.parent_path() / "tmpMsplit";

    createDir(tmpIn);
    createDir(tmpSplit);
    createDir(tmpMsplit);

    //fs::path tmpInPath = tmpIn / fs::path(matched).filename();
    //std::string tmpInNoExt = tmpInPath.string().substr(0,tmpInPath.size()-4);

    std::string tmpInPath = (tmpIn / "tmp").string();

    std::string call = "awk '/^@/ {hdr = hdr $0 ORS; next;}";
    call += "( (++numLines) % " +  std::to_string(entries) + " ) == 1 { if($0 == prev ) {--numLines}";
    call += "else {close(out); out = \"" + tmpInPath + "_\" (++numBlocks) \".sam\"; printf \"%s\", hdr > out; numLines = 1;}}";
    call += "{print > out; prev = $0}' " + matched;

    system(call.c_str());

    // list all files in tmpIn directory
    std::vector<fs::path> subfilesInput; //
    copy(fs::directory_iterator(tmpIn), fs::directory_iterator(), back_inserter(subfilesInput));
    sort(subfilesInput.begin(),subfilesInput.end()); // sort the content

    fs::path splitsPath;
    fs::path msplitsPath;
    std::vector<fs::path> subfilesSplits;
    std::vector<fs::path> subfilesMsplits;

    int counter = 1;
    for(auto &file : subfilesInput) {
        //     std::cout << tmpSplit / file.filename() << std::endl;
//        fs::path splitsPath = tmpSplit / file.stem();
        //       fs::path msplitsPath = tmpMsplit / file.stem();
        subfilesSplits.push_back(tmpSplit / file.filename());
        subfilesMsplits.push_back(tmpMsplit / file.filename());
//        subfilesSplits.push_back(fs::path(splitsPath.string().substr(0,splitsPath.size()-2)+"_splits_"+std::to_string(counter)+".sam"));
        //       subfilesMsplits.push_back(fs::path(msplitsPath.string().substr(0,msplitsPath.size()-2)+"_msplits_"+std::to_string(counter)+".sam"));
        counter++;
    }

    std::vector<std::vector<fs::path>> subfiles;
    subfiles.push_back(subfilesInput);
    subfiles.push_back(subfilesSplits);
    subfiles.push_back(subfilesMsplits);

    return subfiles;
}


// start
void SplitReadCalling::start(pt::ptree sample) {
    // reset counts
    readscount = 0;
    alignedcount = 0;
    splitscount = 0;
    msplitscount = 0;
    nsurvivedcount = 0;

    // input
    pt::ptree input = sample.get_child("input");
    std::string matched = input.get<std::string>("matched");

    int threads = params["threads"].as<int>();
    int entries = std::round(countSamEntries(matched,"") / threads) + threads ;

    // output
    pt::ptree output = sample.get_child("output");
    std::string splits = output.get<std::string>("splits");
    std::string multsplits = output.get<std::string>("multsplits");

    //
    std::vector<std::vector<fs::path>> subfiles = splitInputFile(matched, splits, entries);

# pragma omp parallel for num_threads(params["threads"].as<int>())
    for(unsigned i=0;i<subfiles[0].size();++i) {
        iterate(subfiles[0][i].string(), subfiles[1][i].string(), subfiles[2][i].string());
    }

    // merge results
    std::string callSplits;
    std::string callMsplits;
    if(subfiles[0].size() > 0) {
        for(unsigned i=0;i<subfiles[0].size();++i) {
            if(i==0) {
                callSplits = "awk '/^@/ {print}' ";
                callSplits += subfiles[1][i].string() + " > " + splits;
                callMsplits = "awk '/^@/ {print}' ";
                callMsplits += subfiles[2][i].string() + " > " + multsplits;
                system(callSplits.c_str());
                system(callMsplits.c_str());
            }
            callSplits = "awk '!/^@/ {print}' ";
            callSplits += subfiles[1][i].string() + " >> " + splits;
            callMsplits = "awk '!/^@/ {print}' ";
            callMsplits += subfiles[2][i].string() + " >> " + multsplits;
            system(callSplits.c_str());
            system(callMsplits.c_str());
        }
        // remove tmp folders
        fs::remove_all(subfiles[0][0].parent_path());
        fs::remove_all(subfiles[1][0].parent_path());
        fs::remove_all(subfiles[2][0].parent_path());
    }

    // generate stats
    if(params["stats"].as<std::bitset<1>>() == 1) {
        fs::path statsfile = fs::path(params["outdir"].as<std::string>()) / "detectStat.txt";
        std::fstream fs;
        fs.open(statsfile.string(), std::fstream::app);
        fs << fs::path(matched).stem().string() << "\t";
        fs << countSamEntries(splits,"")/2 << "\t";
        fs << countSamEntries(multsplits," | cut -f1 | sort | uniq ")<< "\t";
        fs << "\n";
        fs.close();
    }
}




/*
SplitReadCalling::SplitReadCalling(po::variables_map params) : params(params), features(params) {
    //
}

SplitReadCalling::~SplitReadCalling() {}


void SplitReadCalling::start(pt::ptree sample) {
    // input
    pt::ptree input = sample.get_child("input");
    std::string matched = input.get<std::string>("matched");

    // get other usefule variables
    int threads = params["threads"].as<int>();

    // output
    pt::ptree output = sample.get_child("output");
    std::string splits = output.get<std::string>("splits");
    std::string multsplits = output.get<std::string>("multsplits");

    // reset counts
    readscount = 0;
    alignedcount = 0;
    splitscount = 0;
    msplitscount = 0;
    nsurvivedcount = 0;


}*/






