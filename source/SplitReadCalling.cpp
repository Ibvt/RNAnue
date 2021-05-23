#include "SplitReadCalling.hpp"

SplitReadCalling::SplitReadCalling(po::variables_map params) : 
    params(params) {
        readCount = 0;
        splitCount = 0;
        multCount = 0;
}

//
void SplitReadCalling::iterate(std::string matched, std::string splits, std::string multsplits) {
    // input .sam record
    seqan3::alignment_file_input fin{matched, 
		seqan3::fields<seqan3::field::id,
						seqan3::field::flag,
                        seqan3::field::ref_id,
						seqan3::field::ref_offset,
						seqan3::field::mapq,
						seqan3::field::cigar,
						seqan3::field::seq,
						seqan3::field::tags,
						seqan3::field::alignment>{}};
    

    std::vector<size_t> ref_lengths{};
	for(auto &info : fin.header().ref_id_info) {
		ref_lengths.push_back(std::get<0>(info));
	}
    std::deque<std::string> ref_ids = fin.header().ref_ids();
   
    // output file for splits
    seqan3::alignment_file_output splitsfile{splits, ref_ids, ref_lengths,
        seqan3::fields<
            seqan3::field::id,
            seqan3::field::flag,
            seqan3::field::ref_id,
            seqan3::field::ref_offset,
            seqan3::field::mapq,
            seqan3::field::cigar,
            seqan3::field::seq,
            seqan3::field::tags>{}};
   
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
	std::list<record_type> splitrecords{};
	
	for (auto & rec : fin) {
//        seqan3::debug_stream << rec << std::endl;
        // head to next read if unmapped
		if(static_cast<bool>(seqan3::get<seqan3::field::flag>(rec) & seqan3::sam_flag::unmapped)) {
			continue;
		}
        if(seqan3::get<seqan3::field::mapq>(rec) <= params["mapquality"].as<int>()) {
            continue;
        }
		
		std::string QNAME = seqan3::get<seqan3::field::id>(rec);
		if((currentQNAME != "") && (currentQNAME != QNAME)) {
            readCount++;
			process(splitrecords, splitsfile, multsplitsfile);	
			splitrecords.clear();
			splitrecords.push_back(rec);
			currentQNAME = QNAME;
		} else {
			splitrecords.push_back(rec);
			currentQNAME = QNAME;
		}
	}
	if(!splitrecords.empty()) {
		process(splitrecords, splitsfile, multsplitsfile);
		splitrecords.clear();
        readCount++;
	}
}


void SplitReadCalling::process(auto &splitrecords, auto &splitsfile, auto &multsplitsfile) {
    Splts splits;

	seqan3::cigar_op cigarOp; // operation of cigar string
	uint32_t cigarOpSize;

    // absolute position in read/alignment (e.g., 1 to end) 
	uint32_t startPosRead; uint32_t endPosRead; 
    // position in split read (e.g., XX:i to XY:i / 14 to 20)
    uint32_t startPosSplit; uint32_t endPosSplit;


//    std::vector<seqan3::cigar> cigarsplit;
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
        auto xctag = tags.get<"XC"_tag>();

        if(xjtag >= 2) { // splits in SAMfile need to consists of at least 2 segments
            seqan3::debug_stream << "------------ original read" << std::endl;
            seqan3::debug_stream << "QNAME: " << qname << std::endl;
            seqan3::debug_stream << "FLAG: " << flag << std::endl;
            seqan3::debug_stream << "CIGAR: " << cigar << std::endl;
            seqan3::debug_stream << "xxtag: " << xxtag << std::endl;
            seqan3::debug_stream << "xytag: " << xytag << std::endl;
            seqan3::debug_stream << "------------" << std::endl;

            // intialise start & end of reads
            startPosRead = 1;
            endPosRead = 0;
            
            // determine the start position within split
            startPosSplit = xxtag;
            endPosSplit = xxtag-1;

            std::cout << "------------- before reading next read -----------" << std::endl;
            std::cout << "splitID: " << splitID << std::endl;
            std::cout << "segmentNR = " << segmentNr << " / total segments =" << xjtag << std::endl;
            std::cout << "--------------------------------------------------" << std::endl;

            // check the cigar string for splits
            for(auto &cig : cigar) {
                // determine size and operator of cigar element
                cigarOpSize = get<uint32_t>(cig);
                cigarOp = get<seqan3::cigar_op>(cig);

                //seqan3::debug_stream << "CIGAR element: " << cigarOpSize << " " << cigarOp << std::endl;
               
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
                    
                    std::cout << "------------- include N -----------" << std::endl;
                    std::cout << "splitID: " << splitID << std::endl;
                    std::cout << "segmentNR = " << segmentNr << " / total segments =" << xjtag << std::endl;
                    std::cout << "--------------------------------------------------" << std::endl;

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
                    seqan3::debug_stream << "\tstartPosRead: " << startPosRead << std::endl;
                    seqan3::debug_stream << "\tendPosRead: " << endPosRead << std::endl;
                }
            }

            seqan3::debug_stream << "seq: " << seq << std::endl;
            auto subSeq = seq | seqan3::views::slice(startPosRead-1,endPosRead);
            seqan3::debug_stream << "subSeq: " << subSeq << std::endl;
                    
            seqan3::sam_tag_dictionary ntags;
            ntags.get<"XX"_tag>() = startPosSplit;
            ntags.get<"XY"_tag>() = endPosSplit;
            ntags["XN"_tag] = splitID;
				
            filterSegments(*it, refOffset, cigarSplit, subSeq, ntags, curated);

            
            segmentNr++;
            
            std::cout << "------------- at the end -----------" << std::endl;
            std::cout << "splitID: " << splitID << std::endl;
            std::cout << "segmentNR = " << segmentNr << " / total segments =" << xjtag << std::endl;
            std::cout << "--------------------------------------------------" << std::endl;

            if(segmentNr == xjtag) {
                std::map<int, std::vector<SamRecord>>::iterator itPutative = putative.begin();
                putative.insert(itPutative, std::make_pair(splitID,curated));
                //curated.clear();
                segmentNr = 0;
                ++splitID;
            }
            ++it;
        } else {
            ++it;
        }
    }

    if(putative.size() > 0) { // split reads found
        std::cout << "putative " <<  putative.size() << std::endl;
        std::vector<SamRecord> p1;
        std::vector<SamRecord> p2;


        // value for complementarity and hybridization energy
        std::pair<float,float> filters;

        std::string p1Qname;
        std::string p2Qname;

        std::map<int, std::vector<SamRecord>>::iterator itSplits = putative.begin();
        for(itSplits;itSplits != putative.end(); ++itSplits) {
            std::vector<SamRecord> splits = itSplits->second;
            if(splits.size() > 1) { // splits consists at least of 2 segments
                for(unsigned i=0;i<splits.size();++i) {
                    for(unsigned j=i+1;j<splits.size();++j) {
                        if(p1.empty() && p2.empty()) {
                            double res = complementarity(seqan3::get<seqan3::field::seq>(splits[i]), seqan3::get<seqan3::field::seq>(splits[j]));
                            float initHyb = hybridize(seqan3::get<seqan3::field::seq>(splits[i]), seqan3::get<seqan3::field::seq>(splits[j])); 
                            filters = std::make_pair(0.6, initHyb);
                            addFilterToSamRecord(splits[i], filters);
                            addFilterToSamRecord(splits[j], filters);
                            splitSegments.push_back(std::make_pair(splits[i],splits[j]));
                        } else {
                            float cmpl = 0.7;
                            float hyb = hybridize(seqan3::get<seqan3::field::seq>(splits[i]), seqan3::get<seqan3::field::seq>(splits[j])); 
                            if(cmpl > filters.first) { // complementarity is higher
                                splitSegments.clear();
                                filters = std::make_pair(cmpl,hyb);
                                addFilterToSamRecord(splits[i], filters);
                                addFilterToSamRecord(splits[j], filters);
                                splitSegments.push_back(std::make_pair(splits[i],splits[j]));
                            } else {
                                if(cmpl == filters.first) { // complementarity is identical
                                    if(hyb > filters.second) { // ... instead check for hybridization
                                        splitSegments.clear();
                                        filters = std::make_pair(cmpl,hyb);
                                        addFilterToSamRecord(splits[i], filters);
                                        addFilterToSamRecord(splits[j], filters);
                                        splitSegments.push_back(std::make_pair(splits[i],splits[j]));
                                    } else {
                                        if(hyb == filters.second) {
                                            addFilterToSamRecord(splits[i], filters);
                                            addFilterToSamRecord(splits[j], filters);
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


    // write to file
    if(splitSegments.size() == 1) {
        writeSamFile(splitsfile, splitSegments);
        splitCount++;
    } else {
        writeSamFile(multsplitsfile, splitSegments);
        multCount++;
    }
}

// write splits back to file
void SplitReadCalling::writeSamFile(auto &samfile, std::vector<std::pair<SamRecord,SamRecord>> splits) {
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

    // filter for mapping quality
    if(seq.size() <= params["minfraglen"].as<int>()) {
    //    return;
    }

    std::cout << "adding" << std::endl;

    // std::cout << "length " << seq.size() << std::endl;
    curated.push_back(segment);
}


void SplitReadCalling::addFilterToSamRecord(SamRecord &rec, std::pair<float,float> filters) {
    seqan3::get<seqan3::field::tags>(rec)["XC"_tag] = filters.first;
    seqan3::get<seqan3::field::tags>(rec)["XE"_tag] = filters.second;
}


double SplitReadCalling::complementarity2(std::span<seqan3::dna5> &seq1, std::span<seqan3::dna5> &seq2) {
    std::cout << "complementarity" << std::endl;
    int m = seq1.size();
    int n = seq2.size();

    int matrix[m+1][n=1];
    // initialise matrix
    for(unsigned i = 0;i<m+1;++i) {
        for(unsigned j = 0;j<n+1;++j) {
            matrix[i][j] = 0; 
        }
    }

    std::vector<seqan3::dna5> seq1_vec;
    std::vector<seqan3::dna5> seq2_vec;

    for(auto &nt : seq1) {
        seq1_vec.push_back(nt);
    }
    for(auto &nt : seq2) {
        seq2_vec.push_back(nt);
    }

    int score = 0;
    std::vector<float> nextEntry{};

    //
    for(unsigned i=1;i<m+1;++i) {
        for(unsigned j=1;j<n+1;++j) {
            //
            
            if((seq1_vec[i-1] == 'A'_dna5 && seq2_vec[j-1] == 'T'_dna5)
                || (seq1_vec[i-1] == 'T'_dna5 && seq2_vec[j-1] == 'A'_dna5)
                || (seq1_vec[i-1] == 'G'_dna5 && seq2_vec[j-1] == 'C'_dna5)
                || (seq1_vec[i-1] == 'C'_dna5 && seq2_vec[j-1] == 'G'_dna5)
                || (seq1_vec[i-1] == 'G'_dna5 && seq2_vec[j-1] == 'U'_dna5)
                || (seq1_vec[i-1] == 'U'_dna5 && seq2_vec[j-1] == 'G'_dna5)) {
                score = 1;
            } else {
                score = -1;
            }

            nextEntry.push_back(0);
            nextEntry.push_back(matrix[i-1][j-1] + score);
            nextEntry.push_back(matrix[i-1][j] - 3);
            nextEntry.push_back(matrix[i][j-1] - 3);

            matrix[i][j] = *max_element(nextEntry.begin(),nextEntry.end());

        }
    }


    // traceback


    float bla;
    return bla;
}


//
double SplitReadCalling::complementarity(std::span<seqan3::dna5> &seq1, std::span<seqan3::dna5> &seq2) {
    std::cout << "start complementarity" << std::endl;
  
    /*



    std::vector<seqan3::dna5> seq1_vec;
    std::vector<seqan3::dna5> seq2_vec;

    for(auto &nt : seq1) {
        seq1_vec.push_back(nt);
    }
    for(auto &nt : seq2) {
        seq2_vec.push_back(nt);
    }
	
    // reverse rna1
	seqan3::dna5_vector rna1Rev;
	for(unsigned z=seq1_vec.size();z-- > 0;) {
		rna1Rev.push_back(seq1_vec[z]);
	}
    
    std::pair p{rna1Rev,seq2_vec};

    auto method = seqan3::align_cfg::method_local{};
    //seqan3::align_cfg::scoring_scheme scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}}};
    // apply specific scoring scheme
    //
    //
    //
    */


    /*
    seqan3::nucleotide_scoring_scheme scheme; // hamming distance is default
    scheme.score('A'_dna15, 'T'_dna15) = 1;
    scheme.score('T'_dna15, 'A'_dna15) = 1;
    scheme.score('G'_dna15, 'C'_dna15) = 1;
    scheme.score('C'_dna15, 'G'_dna15) = 1;
    scheme.score('G'_dna15, 'T'_dna15) = 1;
    scheme.score('A'_dna15, 'G'_dna15) = -1;
    scheme.score('C'_dna15, 'T'_dna15) = -1;
    */
    

    /*
    seqan3::align_cfg::scoring_scheme scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}}};
    seqan3::align_cfg::gap_cost_affine gap_costs{seqan3::align_cfg::open_score{-3}, seqan3::align_cfg::extension_score{-2}};
    
    // Configure the output:
    auto output_config = seqan3::align_cfg::output_score{} |
    	seqan3::align_cfg::output_begin_position{} |
    	seqan3::align_cfg::output_end_position{} |
    	seqan3::align_cfg::output_alignment{};
    
    auto config = method | scheme | gap_costs | output_config;

    auto results = seqan3::align_pairwise(p, config);
    auto & res = *results.begin();
//	seqan3::debug_stream << "Score: " << res.score() << '\n';
//	seqan3::debug_stream << "Begin: (" << res.sequence1_begin_position() << "," << res.sequence2_begin_position() << ")\n";
//	seqan3::debug_stream << "End: (" << res.sequence1_end_position() << "," << res.sequence2_end_position() << ")\n";
//
    
    */

	/*
	// reverse rna1
	seqan3::dna5_vector rna1Rev;
	for(unsigned z=rna1.size();z-- > 0;) {
		rna1Rev.push_back(rna1[z]);
	}

    std::pair p{rna1Rev,rna2};

    // Configure the output:
    auto output_config = seqan3::align_cfg::output_score{} |
    	seqan3::align_cfg::output_begin_position{} |
    	seqan3::align_cfg::output_end_position{} |
    	seqan3::align_cfg::output_alignment{};

    seqan3::nucleotide_scoring_scheme scheme; // hamming distance is default

    
   // scheme.score('A'_dna15, 'T'_dna15) = 1;
    //scheme.score('G'_dna15, 'C'_dna15) = 1;
	//scheme.score('G'_dna15, 'T'_dna15) = 1;
    //scheme.score('A'_dna15, 'G'_dna15) = -1;
    //scheme.score('C'_dna15, 'T'_dna15) = -1;

    // Configure the alignment kernel.
    auto config = seqan3::align_cfg::method_local{} | 
		seqan3::align_cfg::scoring_scheme{scheme} | output_config;

    seqan3::debug_stream << "rna1: " << rna1Rev << std::endl;
    seqan3::debug_stream << "rna2: " << rna2 << std::endl;

    auto results = seqan3::align_pairwise(p, config);
    auto & res = *results.begin();
	seqan3::debug_stream << "Score: " << res.score() << '\n';
	seqan3::debug_stream << "Begin: (" << res.sequence1_begin_position() << "," << res.sequence2_begin_position() << ")\n";
	seqan3::debug_stream << "End: (" << res.sequence1_end_position() << "," << res.sequence2_end_position() << ")\n";

	*/
    
    return 0.0;
}


double SplitReadCalling::hybridize(std::span<seqan3::dna5> &seq1, std::span<seqan3::dna5> &seq2) {
    std::cout << "hyrbidization energy" << std::endl;
	std::string rna1str = "";
	std::string rna2str = "";

	for(unsigned i=0;i<seq1.size();++i) {
		rna1str += seq1[i].to_char();
	}
	for(unsigned i=0;i<seq2.size();++i) {
		rna1str += seq2[i].to_char();
	}

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



// start
void SplitReadCalling::start(pt::ptree sample) {
    // reset 
    readCount = 0;
    splitCount = 0;
    multCount = 0;

    pt::ptree input = sample.get_child("input");
    std::string matched = input.get<std::string>("matched");
    
    pt::ptree output = sample.get_child("output");
    std::string splits = output.get<std::string>("splits");

    std::string multsplits = output.get<std::string>("multsplits");

//    std::cout << matched << std::endl;
 //   std::cout << splits << std::endl;

    // iterate through the reads
    iterate(matched, splits, multsplits);

    // generate stats
    if(params["stats"].as<std::bitset<1>>() == 1) {
        fs::path statsfile = fs::path(params["outdir"].as<std::string>()) / "detectStat.txt";
        std::fstream fs;
        fs.open(statsfile.string(), std::fstream::app);

        fs << fs::path(matched).stem().string() << "\t";
        fs << readCount << "\t";
        fs << splitCount << "\t";
        fs << multCount << std::endl;
        fs.close();
    }
}

