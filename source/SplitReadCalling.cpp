#include "SplitReadCalling.hpp"

SplitReadCalling::SplitReadCalling(po::variables_map params) : 
    params(params) {
}

//
void SplitReadCalling::iterate(std::string matched, std::string splits) {
    // input .sam record
    seqan3::alignment_file_input fin{matched, 
		seqan3::fields<seqan3::field::id,
						seqan3::field::seq,
						seqan3::field::flag,
                        seqan3::field::ref_id,
						seqan3::field::cigar,
						seqan3::field::ref_offset,
						seqan3::field::tags,
						seqan3::field::alignment>{}};
    

    std::vector<size_t> ref_lengths{};
	for(auto &info : fin.header().ref_id_info) {
		ref_lengths.push_back(std::get<0>(info));
	}
    std::deque<std::string> ref_ids = fin.header().ref_ids();
    
    seqan3::alignment_file_output splitsfile{splits, ref_ids, ref_lengths,
        seqan3::fields<
            seqan3::field::id,
            seqan3::field::flag,
            seqan3::field::ref_id,
            seqan3::field::ref_offset,
            seqan3::field::cigar,
            seqan3::field::seq,
            seqan3::field::tags>{}};


    std::string currentQNAME = "";
	//std::vector<std::string> splitrecord; // all split combinations
    using record_type = typename decltype(fin)::record_type;
	std::vector<record_type> splitrecords{};
	
	for (auto & rec : fin) {
//        seqan3::debug_stream << rec << std::endl;
        // head to next read if unmapped
		if(static_cast<bool>(seqan3::get<seqan3::field::flag>(rec) & seqan3::sam_flag::unmapped)) {
			continue;
		}
		//seqan3::debug_stream << "flag: " << seqan3::get<seqan3::field::flag>(rec) << std::endl;
		
		std::string QNAME = seqan3::get<seqan3::field::id>(rec);
		if((currentQNAME != "") && (currentQNAME != QNAME)) {
			process(splitrecords, splitsfile);	
			splitrecords.clear();
			splitrecords.push_back(rec);
			currentQNAME = QNAME;

		} else {
			splitrecords.push_back(rec);
			currentQNAME = QNAME;
		}
	}

	if(!splitrecords.empty()) {
		process(splitrecords, splitsfile);
		splitrecords.clear();
	}
}


void SplitReadCalling::process(auto &splitrecords, auto &splitsfile) {
    Splts splits;

	seqan3::cigar_op cigarOp; // operation of cigar string
	uint32_t cigarOpSize;

	uint32_t startPosRead; 
	uint32_t endPosRead;

    uint32_t startPosSplit;
    uint32_t endPosSplit;

    std::vector<seqan3::cigar> cigSubstr;
    seqan3::sam_tag_dictionary tags;

    // do something with the length
   
    // iterate through each 
    for(unsigned i=0;i<splitrecords.size();++i) {
		std::string qname = "";
		
        // extract the tags (information about split reads) 
        tags = seqan3::get<seqan3::field::tags>(splitrecords[i]);
		auto xhtag = tags.get<"XH"_tag>();
		auto xjtag = tags.get<"XJ"_tag>();
		auto xxtag = tags.get<"XX"_tag>();
		auto xytag = tags.get<"XY"_tag>();

		if(xjtag != 0) { // there is a split in SAM record - number of splits != 0
			qname = seqan3::get<seqan3::field::id>(splitrecords[i]);
            seqan3::sam_flag flag = seqan3::get<seqan3::field::flag>(splitrecords[i]);
			std::optional<int32_t> refID = seqan3::get<seqan3::field::ref_id>(splitrecords[i]);
			std::optional<int32_t> refOffset = seqan3::get<seqan3::field::ref_offset>(splitrecords[i]);
            std::vector<seqan3::cigar> cigar{seqan3::get<seqan3::field::cigar>(splitrecords[i])};
            std::vector<seqan3::cigar> cigarSplit{};
            
            // initialize start and endposition within read
			startPosRead = xxtag; // start equals start in samtag
			endPosRead = xytag;  // end will counted via cigar
           
            // 
            startPosSplit = 1;
            endPosSplit = 0;
		
            // sequence 
            seqan3::dna5_vector seq = seqan3::get<seqan3::field::seq>(splitrecords[i]);

            //seqan3::debug_stream << seq[0] << std::endl;
            uint32_t seqLen = seq.size();


            /*
            seqan3::debug_stream << qname << std::endl;
            seqan3::debug_stream << "complete sequence: " << (seq | seqan3::views::to_char) << std::endl;
            seqan3::debug_stream << "sequence length: " << seqLen << std::endl;

            seqan3::debug_stream << "xxtag: " << xxtag << " ";
            seqan3::debug_stream << "xytag: " << xytag << std::endl;
            seqan3::debug_stream << "cigar string: " << cigar << std::endl;
            seqan3::debug_stream << "go through cigar" << std::endl;
            */
            
            int basesPassed = 0;
            
            for(auto &cig : cigar) {
                //voegel
                cigarOpSize = get<uint32_t>(cig);
                cigarOp = get<seqan3::cigar_op>(cig);
                std::cout << cigarOpSize << seqan3::to_char(cigarOp);
               
                if(cigarOp != 'N'_cigar_op) {
                    seqan3::cigar cigEl{cigarOpSize, cigarOp};
                    cigarSplit.push_back(cigEl);
                }

                if(cigarOp == 'N'_cigar_op) {
                    seqan3::debug_stream << "includes N" << std::endl;
                    seqan3::debug_stream << "startPosSplit: " << startPosSplit << " endPosSplit: " << endPosSplit << std::endl;

                    // determine sequence of the splits
                    // auto splitSeq = seq | seqan3::views::slice(startPosSplit-1,endPosSplit);
                    //auto splitSeq = seq | seqan3::views::slice(startPosRead-1,endPosRead);

                    // 
                    seqan3::dna5_vector splitSeq;
                    for(uint32_t z=startPosSplit-1;z<endPosSplit;z++) {
                        splitSeq.push_back(seq[z]);
                    }
                    seqan3::debug_stream << "subsequence: " << splitSeq << std::endl;
                   
                    seqan3::sam_tag_dictionary dict{};
                    dict.get<"XX"_tag>() = startPosSplit;
                    dict.get<"XY"_tag>() = endPosSplit;

                    seqan3::debug_stream << "cigar of split" << cigarSplit << std::endl;

                    splits.push_back(
                            std::make_tuple(
                                qname, flag, refID, refOffset, cigarSplit, splitSeq, dict));

                    cigarSplit.clear();
                   

					startPosSplit = endPosSplit+1;

//                    basesPassed += cigarOpSize;
 //                   refOffset.emplace(refOffset.value()+basesPassed);
                    //refOffset += endPosRead;

                } else {
                    basesPassed += cigarOpSize;
                    if(cigarOp != 'D'_cigar_op) {
                        endPosSplit += cigarOpSize;
                    }
                }
            }
            seqan3::debug_stream << "\nend of read reached: start: " << startPosRead + startPosSplit -1 << " end: " << startPosRead + endPosSplit -1 << std::endl; 
//            auto splitSeq = seq | seqan3::views::slice(startPosRead,endPosRead);
            
            seqan3::dna5_vector splitSeq2;
            for(uint32_t z=startPosSplit-1;z<endPosSplit;z++) {
                splitSeq2.push_back(seq[z]);
            }
            seqan3::debug_stream << "split sequence: " << splitSeq2 << std::endl;
                   
			
            seqan3::sam_tag_dictionary dict{};
            dict.get<"XX"_tag>() = startPosRead;
            dict.get<"XY"_tag>() = endPosRead;

            splits.push_back(
                std::make_tuple(
                    qname, flag, refID, refOffset, cigarSplit, splitSeq2, dict));

        }
    }

    if(splits.size() > 1) {
        std::cout << "\n\n############### output" << std::endl;
        std::cout << splits.size() << std::endl;

		
        for(unsigned c=0;c<splits.size();++c) {
            for(unsigned d=c+1;d<splits.size();++d) {

				double cmpl = complementarity(std::get<5>(splits[c]),std::get<5>(splits[d]));
				std::cout << "complementarity: " << cmpl << std::endl;
				std::cout << "\n";

				double nrg = hybridize (std::get<5>(splits[c]),std::get<5>(splits[d]));
				
				seqan3::sam_tag_dictionary dict2nd = std::get<6>(splits[d]);
                dict2nd.get<"FC"_tag>() = std::to_string(cmpl);
                dict2nd.get<"FE"_tag>() = std::to_string(nrg);
//				dict2nd.get<"AB"_tag>() = cmpl;
//				dict2nd.get<"DD"_tag>() = nrg;

				splitsfile.emplace_back(
					std::get<0>(splits[c]),
					std::get<1>(splits[c]),
					std::get<2>(splits[c]),
					std::get<3>(splits[c]),
					std::get<4>(splits[c]),
					std::get<5>(splits[c]),
					dict2nd);

				splitsfile.emplace_back(
					std::get<0>(splits[d]),
					std::get<1>(splits[d]),
					std::get<2>(splits[d]),
					std::get<3>(splits[d]),
					std::get<4>(splits[d]),
					std::get<5>(splits[d]),
					std::get<6>(splits[d]));
                    
            }
        } 
    }
}

//
double SplitReadCalling::complementarity(seqan3::dna5_vector rna1, seqan3::dna5_vector rna2) {
    std::cout << "start complementarity" << std::endl;


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
    
    return 0.0;
}


double SplitReadCalling::hybridize(seqan3::dna5_vector rna1, seqan3::dna5_vector rna2) {
	std::cout << "hyrbidization energy" << std::endl;

	std::string rna1str = "";
	std::string rna2str = "";

	for(unsigned i=0;i<rna1.size();++i) {
		rna1str += rna1[i].to_char();
	}
	for(unsigned i=0;i<rna2.size();++i) {
		rna1str += rna2[i].to_char();
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
    pt::ptree input = sample.get_child("input");
    std::string matched = input.get<std::string>("matched");
    
    pt::ptree output = sample.get_child("output");
    std::string splits = output.get<std::string>("splits");

    std::cout << matched << std::endl;
    std::cout << splits << std::endl;

    // iterate through the reads
    iterate(matched, splits);



}


