#include "Align.hpp"


Align::Align(po::variables_map params) : 
    params(params) {

    std::cout << "create align object" << std::endl;

    buildIndex();
}

void Align::alignReads(std::string query, std::string matched, std::string unmatched) {
    std::string align = "segemehl.x";
    align += " -S ";
    align += " -A " + std::to_string(params["accuracy"].as<int>()); 
    align += " -U " + std::to_string(params["minfragsco"].as<int>());
    align += " -W " + std::to_string(params["minsplicecov"].as<int>());
    align += " -Z " + std::to_string(params["minfraglen"].as<int>());
    align += " -t " + std::to_string(params["threads"].as<int>());
    align += " -i " + index;
    align += " -d " + params["dbref"].as<std::string>();
    align += " -q " + query;
    align += " -o " + matched;

    const char* call = align.c_str();
    system(call);

    std::cout << align << std::endl;
}

void Align::buildIndex() {
    // retrieve path of reference genome
    std::string ref = params["dbref"].as<std::string>();

    fs::path outDir = fs::path(params["outdir"].as<std::string>());
    fs::path gen = outDir / fs::path(ref).replace_extension(".idx").filename();

    if(fs::exists(gen)) {
        std::cout << "segemehl index found on filesystem\n";
    } else {
        std::cout << "generate index " << "\n";
        std::string genIndex = "segemehl.x -x " + gen.string() + " -d " + ref;
        const char* call = genIndex.c_str();
        system(call);
    }
    index = gen.string();
}

void Align::detSplits(std::string matched, std::string splits) {
    std::cout << "determine the split reads\n";

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
	
    //seqan3::debug_stream << fin.header().ref_dict << '\n';
	//std::ofstream splitsFile(splits);
    seqan3::alignment_file_output splitsfile{splits, ref_ids, ref_lengths,
        seqan3::fields<
			seqan3::field::id,
            seqan3::field::flag,
			seqan3::field::ref_id,
			seqan3::field::ref_offset,
            seqan3::field::cigar,
            seqan3::field::seq,
            seqan3::field::tags>{}};

	// copy dictionary with ref_ids	
	//splitsfile.header().ref_dict = refIDs;

    std::string currentQNAME = "";
	//std::vector<std::string> splitrecord; // all split combinations
    using record_type = typename decltype(fin)::record_type;
	std::vector<record_type> splitrecords{};
	
	for (auto & rec : fin) {
		if(static_cast<bool>(seqan3::get<seqan3::field::flag>(rec) & seqan3::sam_flag::unmapped)) {
			continue;
		}
		//seqan3::debug_stream << "flag: " << seqan3::get<seqan3::field::flag>(rec) << std::endl;
		
		std::string QNAME = seqan3::get<seqan3::field::id>(rec);
		if((currentQNAME != "") && (currentQNAME != QNAME)) {
			processSplits(splitrecords, splitsfile);	
			splitrecords.clear();
			splitrecords.push_back(rec);
			currentQNAME = QNAME;

		} else {
			splitrecords.push_back(rec);
			currentQNAME = QNAME;
		}
	}

	if(!splitrecords.empty()) {
		processSplits(splitrecords, splitsfile);
		splitrecords.clear();
	}
}

void Align::processSplits(auto &splitrecords, auto &splitsfile) {
	// QNAME, <leftPosRead,rightPosRead>, 
	Splts splits;

	seqan3::cigar_op cigarOp; // operation of cigar string
	uint32_t cigarOpSize;

	uint32_t startPosRead; 
	uint32_t endPosRead;

    uint32_t startPosSplit;
    uint32_t endPosSplit;

    std::vector<seqan3::cigar> cigSubstr;

    seqan3::sam_tag_dictionary tags;

    for(unsigned i=0;i<splitrecords.size();++i) {
		std::string qname = "";

		tags = seqan3::get<seqan3::field::tags>(splitrecords[i]);
		auto xhtag = tags.get<"XH"_tag>();
		auto xjtag = tags.get<"XJ"_tag>();
		auto xxtag = tags.get<"XX"_tag>();
		auto xytag = tags.get<"XY"_tag>();

		if(xhtag != 0) { // there is a split in SAM record
			qname = seqan3::get<seqan3::field::id>(splitrecords[i]);
            seqan3::sam_flag flag = seqan3::get<seqan3::field::flag>(splitrecords[i]);
			std::optional<int32_t> refID = seqan3::get<seqan3::field::ref_id>(splitrecords[i]);
			std::optional<int32_t> refOffset = seqan3::get<seqan3::field::ref_offset>(splitrecords[i]);
            std::vector<seqan3::cigar> cigar{seqan3::get<seqan3::field::cigar>(splitrecords[i])};
            std::vector<seqan3::cigar> cigarSplit{};

            // initialize start and endposition within read
			startPosRead = xxtag; // start equals start in samtag
			endPosRead = xxtag-1;  // end will counted via cigar

            startPosSplit = 1;
            endPosSplit = 0;

			seqan3::dna5_vector seq = seqan3::get<seqan3::field::seq>(splitrecords[i]);
            uint32_t seqLen = seq.size();

            seqan3::debug_stream << qname << std::endl;
            //seqan3::debug_stream << flag << std::endl;
            seqan3::debug_stream << "complete sequence: " << (seq | seqan3::views::to_char) << std::endl;
            seqan3::debug_stream << "sequence length" << seqLen << std::endl;

            seqan3::debug_stream << "xxtag: " << xxtag << " ";
            seqan3::debug_stream << "xytag: " << xytag << std::endl;
            seqan3::debug_stream << cigar << std::endl;

			for(auto &cig: cigar) {
				cigarOpSize = get<uint32_t>(cig);
				cigarOp = get<seqan3::cigar_op>(cig);
//                std::cout << cigarOpSize << seqan3::to_char(cigarOp);
				if(cigarOp == 'N'_cigar_op) {
                    seqan3::debug_stream << "includes N" << std::endl;
                    seqan3::debug_stream << "startPosRead: " << startPosRead << " endPosRead: " << endPosRead << std::endl;

                    // determine sequence of split
                    seqan3::debug_stream << "determine sequence of split: start: " << startPosSplit  << " end: " << endPosSplit << std::endl;
                    std::span<seqan3::dna5, -1> splitSeq = seq | seqan3::views::slice(startPosSplit-1,endPosSplit);
                    seqan3::debug_stream << "sequence of split: " << splitSeq << std::endl;
                    
                    
                   // tags.get<"XY"_tag>() = endPosRead; 

					splits.push_back(
                            std::make_tuple(
                                qname, 
                                flag, 
                                refID, 
                                refOffset, 
                                cigarSplit, 
                                spanToVec(splitSeq), 
                                tags));

                    cigarSplit.clear();


					startPosRead = ++endPosRead;
					endPosRead = startPosRead-1;

                    // next split
                    startPosSplit = ++endPosSplit; // start
                    endPosSplit = startPosSplit-1; // end
                    tags.get<"XX"_tag>() = startPosRead; // start within read
                    refOffset = *refOffset + cigarOpSize + startPosRead;
                    seq = spanToVec(seq | seqan3::views::slice(endPosSplit, seqLen));
                  

                    seqan3::debug_stream << "rest of sequence : " << seq << std::endl;
                    seqan3::debug_stream << "endPosRead: " << endPosRead << std::endl;
					continue;
				} else {
					if(cigarOp != 'D'_cigar_op) {
                        cigarSplit.push_back(cig);
                        endPosRead += cigarOpSize;
                        endPosSplit += cigarOpSize;

					}
				}
			}
            tags.get<"XY"_tag>() = endPosRead;
//			std::cout << std::endl;
            seqan3::debug_stream << "startPosRead: " << startPosRead;
            seqan3::debug_stream << " endPosRead: " << endPosRead << std::endl;
            seqan3::debug_stream << seq << std::endl;
            splits.push_back(
                    std::make_tuple(
                        qname,
                        flag, 
                        refID,
                        refOffset, 
                        cigarSplit, 
                        seq, 
                        tags));
		}
	}

    
    if(splits.size() >= 2) {
        // 
        std::map<std::pair<uint32_t,uint32_t>,std::pair<double,double>> putSplits;

        seqan3::debug_stream << "##########################" << std::endl;
        std::cout << splits.size() << std::endl;
        for(unsigned i=0;i<splits.size();++i) {
            for(unsigned j=i+1;j<splits.size();++j) {
                auto tagsT1 = std::get<6>(splits[i]);
                auto tagsT2 = std::get<6>(splits[j]);

                auto startT1 = tagsT1.get<"XX"_tag>();
                auto endT1 = tagsT1.get<"XY"_tag>();

                auto startT2 = tagsT2.get<"XX"_tag>();
                auto endT2 = tagsT2.get<"XY"_tag>();

                if((startT1 < startT2) && (endT1 < startT2) || 
                        (startT1 > startT2) && (endT2 < startT1))  {

                    seqan3::debug_stream << "------------new entry--------------" << std::endl;
                    seqan3::debug_stream << "startT1: " << startT1 << std::endl;
                    seqan3::debug_stream << "endT1: " << endT1 << std::endl;
                    seqan3::debug_stream << "startT2: " << startT2<< std::endl;
                    seqan3::debug_stream << "endT2: " << endT2 << std::endl;

                    auto seq1 = (std::get<5>(splits[i]) | seqan3::views::to_char);
                    auto seq2 = (std::get<5>(splits[j]) | seqan3::views::to_char);

                    seqan3::debug_stream << "seq1: " << seq1 << std::endl;
                    seqan3::debug_stream << "seq2: " << seq2 << std::endl;

                    std::cout << std::string(seq1.begin(),seq1.end()) << "&" << std::string(seq2.begin(),seq2.end()) << std::endl;
                    seqan3::debug_stream << "-------------------------------------" << std::endl;

                    hybridize(std::string(seq1.begin(),seq1.end()),std::string(seq2.begin(),seq2.end()));

                    putSplits.insert(
                            std::make_pair(
                                std::make_pair(i,j),
                                std::make_pair(0.0,
                                    hybridize(
                                        std::string(seq1.begin(),seq1.end()),
                                        std::string(seq2.begin(),seq2.end())))));

                   // system(call);
                }
            }
        }
        splits.clear();

        for(auto &spl : putSplits) {
            seqan3::debug_stream << "(" << spl.first.first << "," <<  spl.first.second << ") -> " << spl.second.second << std::endl;


        }


/*

        for(auto &splt : splits) {
            tags = std::get<6>(splt);
            tags.get<"XH"_tag>() = 1; // splits in sam record
            tags.get<"XJ"_tag>() = 2; // splits in read


            splitsfile.emplace_back(
				std::get<0>(splt),
				std::get<1>(splt), 
				std::get<2>(splt),
				std::get<3>(splt),
                std::get<4>(splt),
                std::get<5>(splt),
                tags);
        }*/
    }

}

double Align::hybridize(std::string rna1, std::string rna2) {
    std::string hyb = "echo '" + rna1 + "&" + rna2 + "' | RNAcofold";
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

//
double complementarity(std::string rna1, std::string rna2) {
    std::cout << "start complementarity" << std::endl;


    // perform alignment

    // instatiate matrix
    //


    // fill




    // backtracking
}


//
void Align::start(pt::ptree sample) {
    pt::ptree input = sample.get_child("input");
    pt::ptree output = sample.get_child("output");

    if(params["readtype"].as<std::string>() == "SE") {
        std::string forward = input.get<std::string>("forward");
        std::string matched = output.get<std::string>("matched");
        std::string splits = output.get<std::string>("splits");

        alignReads(forward, matched, splits);
        std::cout << "filter split reads: " << splits << '\n';
        detSplits(matched, splits);
    }
}
        
std::vector<seqan3::dna5> Align::spanToVec(std::span<seqan3::dna5,-1> seq) {
    std::vector<seqan3::dna5> out;
    for(auto &nt : seq) {
        out.push_back(nt);
    }
    return(out);
}

