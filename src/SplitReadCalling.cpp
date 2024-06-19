#include "SplitReadCalling.hpp"

template <typename T>
SafeQueue<T>::SafeQueue() {}

template <typename T>
SafeQueue<T>::~SafeQueue() {}

template <typename T>
void SafeQueue<T>::push(T value) {
    std::lock_guard<std::mutex> lock(mtx);
    q.push(std::move(value));
    cv.notify_one();
}

template <typename T>
bool SafeQueue<T>::pop(T& result) {
    std::unique_lock<std::mutex> lock(mtx);
    cv.wait(lock, [this]() { return !q.empty() || done; });
    if (q.empty()) return false;
    result = std::move(q.front());
    q.pop();
    return true;
}

template <typename T>
void SafeQueue<T>::reset() {
    std::lock_guard<std::mutex> lock(mtx);
    std::queue<T> empty;
    std::swap(q, empty);
    done = false;
}

template <typename T>
void SafeQueue<T>::setDone() {
    std::lock_guard<std::mutex> lock(mtx);
    done = true;
    cv.notify_all();
}

SplitReadCalling::SplitReadCalling(po::variables_map params) {
    this->params = params;
    // create IBPTree (if splicing need to be considered)
    if(params["splicing"].as<std::bitset<1>>() == std::bitset<1>("1")) {
        // order of the IBPTree
        int k = 3; // can be later extracted from config file
        // create IBPT
        this->features = IBPTree(params, k);
    }
    this->condition = ""; // stores the current condition
    this->refIds = std::deque<std::string>(); // stores the reference IDs

    // initialize stats
    this->stats = std::make_shared<Stats>();
    this->replPerCond = -1; // number of replicates per condition

}

SplitReadCalling::~SplitReadCalling() {}

void SplitReadCalling::start(pt::ptree sample, pt::ptree condition) {
    // input
    pt::ptree input = sample.get_child("input");
    std::string matched = input.get<std::string>("matched");

    if(this->condition != "" && this->condition != condition.get<std::string>("condition")) {
        replPerCond = -1; // 'new' condition - reset replicates counter
    }
    this->condition = condition.get<std::string>("condition"); // store current condition (e.g. rpl_37C)
    this->stats->reserveStats(this->condition, ++replPerCond); // increase the number of replicates per condition

    // output
    pt::ptree output = sample.get_child("output");
    std::string single = output.get<std::string>("single");
    std::string splits = output.get<std::string>("splits");
    std::string multsplits = output.get<std::string>("multsplits");

    std::cout << helper::getTime() << "Split Read Calling on sample: " << matched << std::endl;

    srand(time(NULL)); // seed the random number generator
    int randNum = rand() % 1000;
    std::string sortedBam = this->condition + "_" + std::to_string(randNum) + ".bam";
    fs::path sortedBamFile = fs::path(params["outdir"].as<std::string>()) / fs::path("tmp") / fs::path(sortedBam);
    std::string sortedBamFileStr = sortedBamFile.string();
    sort(matched, sortedBamFileStr);

    seqan3::sam_file_input inBam{sortedBamFileStr}; // initialize input file
    std::vector<size_t> refLengths;
    for(auto &info : inBam.header().ref_id_info) {
        refLengths.push_back(std::get<0>(info));
    }
    this->refIds = inBam.header().ref_ids();
    using recordType = typename decltype(inBam)::record_type;

    // create output files
    seqan3::sam_file_output singleOut{single, this->refIds, refLengths}; // initialize output file
    seqan3::sam_file_output splitsOut{splits,  this->refIds, refLengths}; // initialize output file
    seqan3::sam_file_output multsplitsOut{multsplits, this->refIds, refLengths}; // initialize output file
    std::mutex singleOutMutex, splitsOutMutex, multsplitsOutMutex; // create mutex objects for output

    queue.reset(); // reset the queue

    std::vector<std::thread> consumers; // intialize the consumers
    for(int i=0;i<params["threads"].as<int>();++i) {
        consumers.emplace_back([this, &singleOut, &splitsOut, &multsplitsOut,
                                &singleOutMutex, &splitsOutMutex, &multsplitsOutMutex] () {
            this->consumer(singleOut, splitsOut, multsplitsOut, singleOutMutex, splitsOutMutex, multsplitsOutMutex);
        });
    }

    producer(inBam); // start the producer (and add the recordsreads to the queue)

    for(auto& consumer : consumers) {
        consumer.join();
    }
}

void SplitReadCalling::producer(seqan3::sam_file_input<>& inputfile) {
    using record_type = typename seqan3::sam_file_input<>::record_type;
    std::vector<record_type> readRecords;
    std::string QNAME = "", currentQNAME = "";
    for(auto && record : inputfile) {
        // ignore reads if unmapped and do not suffice length requirements
        if(static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped)) { continue; }
        if(static_cast<int>(record.sequence().size()) <= params["minlen"].as<int>()) { continue; }
        if(static_cast<int>(record.mapping_quality()) <= params["mapq"].as<int>()) { continue; }

        // group reads by QNAME
        QNAME = record.id();
        if((currentQNAME != "") && (currentQNAME != QNAME)) {
            this->stats->setReadsCount(this->condition, replPerCond, 1);
            queue.push(readRecords);
            readRecords.clear();
        }
        // prepare for the next iteration
        readRecords.push_back(record);
        currentQNAME = QNAME;
    }
    if(!readRecords.empty()) {
        queue.push(readRecords);
    }
    queue.setDone();
}

void SplitReadCalling::consumer(dtp::BAMOut& singleOut, dtp::BAMOut& splitsOut, dtp::BAMOut& multsplitsOut,
                                std::mutex& singleOutMutex, std::mutex& splitsOutMutex, std::mutex& multsplitsOutMutex) {
    using recordType = typename seqan3::sam_file_input<>::record_type;
    std::vector<recordType> readrecords;
    while(true) {
        if(!queue.pop(readrecords)) {
            break;
        }
        try {
            this->process(readrecords, singleOut, splitsOut, multsplitsOut,
                          singleOutMutex, splitsOutMutex, multsplitsOutMutex);
        } catch (const std::exception& e) {
            std::cerr << helper::getTime() << "Error: " << e.what() << std::endl;
        }
    }
}

template <typename T>
void SplitReadCalling::process(std::vector<T>& readrecords, auto& singleOut,
                               auto& splitsOut, auto& multsplitsOut,
                               auto& singleOutMutex, auto& splitsOutMutex,
                               auto& multsplitsOutMutex) {

    // define variables
    std::string qname = "";
    seqan3::sam_flag flag = seqan3::sam_flag::none;
    std::optional<int32_t> ref_id, ref_offset = std::nullopt;
    std::optional<uint8_t> mapq = std::nullopt;
    dtp::CigarVector cigar = {};
    dtp::DNAVector seq = {};
    seqan3::sam_tag_dictionary tags = {};

    int segNum, splitId = 0; // number of segments and splitID in the read

    // cigar
    uint32_t cigarSize;
    seqan3::cigar::operation cigarOp;
    dtp::CigarVector cigarSplit{}; // cigar for each split
    uint32_t cigarMatch{}; // cigar for the matching part of the split

    uint32_t startPosRead, endPosRead; // absolute position in read/alignment (e.g., 1 to end)
    uint32_t startPosSplit, endPosSplit; // position in split read (e.g., XX:i to XY:i / 14 to 20)

    dtp::DNASpan subSeq;
    dtp::QualSpan subQual;

    std::pair<uint32_t, uint32_t> junction; // splice junction

    // stores the putative splits/segments
    std::vector<SAMrecord> curated;
    std::map<int, std::vector<SAMrecord>> putative{};

    for(auto& rec : readrecords) {
        qname = rec.id();
        flag = rec.flag();
        ref_id = rec.reference_id();
        ref_offset = rec.reference_position();
        mapq = rec.mapping_quality();
        seq = rec.sequence();

        // CIGAR string
        cigar = rec.cigar_sequence();
        cigarSplit.clear();

        tags = rec.tags();
        auto xxtag = tags.get<"XX"_tag>();
        auto xytag = tags.get<"XY"_tag>();
        auto xjtag = tags.get<"XJ"_tag>();
        auto xhtag = tags.get<"XH"_tag>();

        if(tags.get<"XJ"_tag>() < 2 ) { // ignore reads with less than 2 splits (or no splits at all)
            {
                std::lock_guard<std::mutex> lock(singleOutMutex);
                singleOut.push_back(rec);
                this->stats->setAlignedCount(this->condition, replPerCond, 1);
            }
            break;
        }
        startPosRead = 1; // initialize start/end of read
        endPosRead = 0;
        startPosSplit = xxtag;
        endPosSplit = xytag;

        cigarMatch = 0; // reset matches count in cigar
        // check the cigar string for splits
        for(auto& cig: cigar) {
            // determine the size and operator of the cigar element
            cigarSize = get<uint32_t>(cig);
            cigarOp = get<seqan3::cigar::operation>(cig);

            if(cigarOp == 'N'_cigar_operation) {
                dtp::Interval spliceSite = {ref_offset.value() + endPosRead, ref_offset.value() + cigarSize + endPosRead};
                if(!matchSpliceSites(spliceSite, ref_id)) { // if it does not match the splice site
                    // store split
                    seqan3::sam_tag_dictionary newTags{};
                    newTags.get<"XX"_tag>() = startPosSplit;
                    newTags.get<"XY"_tag>() = endPosSplit;
                    newTags.get<"XN"_tag>() = splitId;

                    subSeq = seq | seqan3::views::slice(startPosRead - 1, endPosRead);
                    subQual = rec.base_qualities() | seqan3::views::slice(startPosRead - 1, endPosRead);
                    if(filter(subSeq, cigarMatch)) {
                        storeSegments(rec, ref_offset, subSeq, subQual, cigarSplit, newTags, curated);
                    }

                    // settings to prepare for the next split
                    startPosSplit = endPosSplit+1;
                    startPosRead = endPosRead+1;
                    cigarSplit.clear(); // new split - new CIGAR
                    ref_offset.value() += cigarSize + endPosRead + 1; // adjust left-most mapping position
                    splitId++; // increase split ID
                } else { // matches the splice junction
                    xjtag--; // decrease the number of splits (basically merge the splits)
                }
            } else {
                if(cigarOp == 'S'_cigar_operation &&
                    params["exclclipping"].as<std::bitset<1>>() == std::bitset<1>("1")) {
                    // exclude soft clipping from the alignments
                    if(cigarSplit.size() == 0) { // left soft clipping
                        startPosRead += cigarSize;
                        endPosRead += cigarSize;
                        startPosSplit += cigarSize;
                        endPosSplit += cigarSize;
                    }
                } else {
                    if(cigarOp != 'D'_cigar_operation) {
                        endPosRead += cigarSize;
                        endPosSplit += cigarSize;
                        if(cigarOp == '='_cigar_operation) { cigarMatch += cigarSize; }
                    }
                    cigarSplit.push_back(cig);
                }
            }
        }

        auto subSeq = seq | seqan3::views::slice(startPosRead - 1, endPosRead);
        seqan3::sam_tag_dictionary newTags{};
        newTags.get<"XX"_tag>() = startPosSplit;
        newTags.get<"XY"_tag>() = endPosSplit;
        newTags.get<"XN"_tag>() = splitId;

        subSeq = seq | seqan3::views::slice(startPosRead - 1, endPosRead);
        subQual = rec.base_qualities() | seqan3::views::slice(startPosRead - 1, endPosRead);
        if(filter(subSeq, cigarMatch)) {
            storeSegments(rec, ref_offset, subSeq, subQual, cigarSplit, newTags, curated);
        }
        putative.insert(std::make_pair(splitId, curated));
        curated.clear();
        segNum = 0;
        ++splitId;
    }
    decide(putative, splitsOut, multsplitsOut, splitsOutMutex, multsplitsOutMutex);
}

void SplitReadCalling::decide(std::map<int, std::vector<SAMrecord>>& putative, auto& splitsOut,
                              auto& multsplitsOut, auto& splitsOutMutex,
                              auto& multsplitsOutMutex) {
    std::vector<SAMrecord> p1, p2;
    seqan3::sam_tag_dictionary p1Tags, p2Tags;
    std::pair<double, double> filters;
    std::string p1QNAME, p2QNAME;
    std::vector<std::pair<SAMrecord, SAMrecord>> finalSplits;

    if(putative.size() > 0) {
        for (auto &[splitId, records]: putative) {
            if (records.size() > 1) { // splits contain more than one segment
                for (unsigned i = 0; i < records.size(); ++i) {
                    for (unsigned j = i + 1; j < records.size(); ++j) {
                        p1Tags = records[i].tags();
                        p2Tags = records[j].tags();

                        auto p1Start = p1Tags.get<"XX"_tag>();
                        auto p1End = p1Tags.get<"XY"_tag>();
                        auto p2Start = p2Tags.get<"XX"_tag>();
                        auto p2End = p2Tags.get<"XY"_tag>();

                        // prevent overlap between read position -> same segment of read
                        if ((p1Start >= p1Start && p2Start <= p1End) ||
                            (p1Start >= p2Start && p1Start <= p2End)) {
                            continue;
                        }

                        // determine complementarity/hybridization
                        if (p1.empty() && p2.empty()) {
                            TracebackResult initCmpl = complementarity(records[i].sequence(), records[j].sequence());
                            std::pair<double,std::string> initHyb = hybridization(records[i].sequence(), records[j].sequence());
                            if (!initCmpl.a.empty()) { // split read survived complementarity/sitelenratio cutoff
                                if (initHyb.first <= params["nrgmax"].as<double>()) {
                                    filters = std::make_pair(initCmpl.cmpl, initHyb.first);
                                    addComplementarityToSamRecord(records[i], records[j], initCmpl);
                                    addHybEnergyToSamRecord(records[i], records[j], initHyb);
                                    finalSplits.push_back(std::make_pair(records[i], records[j]));
                                }
                            }
                        } else {
                            TracebackResult cmpl = complementarity(records[i].sequence(), records[j].sequence());
                            std::pair<double,std::string> hyb = hybridization(records[i].sequence(), records[j].sequence());
                            if (!cmpl.a.empty()) { // split read survived complementarity/sitelenratio cutoff
                                if (cmpl.cmpl > filters.first) { // complementarity is higher
                                    finalSplits.clear();
                                    filters = std::make_pair(cmpl.cmpl, hyb.first);
                                    addComplementarityToSamRecord(records[i], records[j], cmpl);
                                    addHybEnergyToSamRecord(records[i], records[j], hyb);
                                    finalSplits.push_back(std::make_pair(records[i], records[j]));
                                } else {
                                    if (cmpl.cmpl == filters.first) {
                                        if (hyb.first > filters.second) {
                                            finalSplits.clear();
                                            filters = std::make_pair(cmpl.cmpl, hyb.first);
                                            addComplementarityToSamRecord(records[i], records[j], cmpl);
                                            addHybEnergyToSamRecord(records[i], records[j], hyb);
                                            finalSplits.push_back(std::make_pair(records[i], records[j]));
                                        }
                                    } else {
                                        if (hyb.first == filters.second) { // same compl & hyb --> ambiguous
                                            addComplementarityToSamRecord(records[i], records[j], cmpl);
                                            addHybEnergyToSamRecord(records[i], records[j], hyb);
                                            finalSplits.push_back(std::make_pair(records[i], records[j]));
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
    if(finalSplits.size() == 1) {
        {
            std::lock_guard<std::mutex> lock(splitsOutMutex);
            writeSAMrecordToBAM(splitsOut, finalSplits);
            this->stats->setSplitsCount(this->condition, replPerCond, 1);
        }
    } else { // multisplits detected
        if(finalSplits.size() > 1) {
            {
                std::lock_guard<std::mutex> lock(multsplitsOutMutex);
                writeSAMrecordToBAM(multsplitsOut, finalSplits);
                this->stats->setMultSplitsCount(this->condition, replPerCond, 1);
            }
        }
    }
}

void SplitReadCalling::writeSAMrecordToBAM(auto& bamfile,
                                           std::vector<std::pair<SAMrecord, SAMrecord>>& records) {
    for(unsigned i=0;i<records.size();++i) {
        bamfile.push_back(records[i].first);
        bamfile.push_back(records[i].second);
    }
}

bool SplitReadCalling::filter(auto& sequence, uint32_t cigarmatch) {
    if(sequence.size() < params["minlen"].as<int>()) { return false; }
    return true;
}

bool SplitReadCalling::matchSpliceSites(dtp::Interval& spliceSites, std::optional<uint32_t> refId) {
    if(params["splicing"].as<std::bitset<1>>() == std::bitset<1>("1")) {
        std::vector<std::pair<Node*,IntervalData*>> ovlps = features.search(this->refIds[refId.value()], spliceSites);
        for(int i=0;i<ovlps.size();++i) { // iterate over all overlapping intervals
            dtp::SpliceJunctions junctions = ovlps[i].second->getJunctions();
            for(auto& [name, junc] : junctions) {
                for(int j=1;j<junc.size()-1;++j) {
                    if(helper::withinRange(spliceSites.first, junc[j-1].second, 5) &&
                    helper::withinRange(spliceSites.second, junc[j].first, 5)) {
                        return true;
                    }
                }
            }
        }
    } else {
        return false;
    }
}

void SplitReadCalling::storeSegments(auto& splitrecord, std::optional<int32_t> refPos, dtp::DNASpan& seq,
                                     dtp::QualSpan& qual, dtp::CigarVector& cigar,
                                     seqan3::sam_tag_dictionary& tags, std::vector<SAMrecord>& curated) {

    SAMrecord segment{}; // create new SAM Record
    segment.id() = splitrecord.id();
    seqan3::sam_flag flag{0};
    if(static_cast<bool>(splitrecord.flag() & seqan3::sam_flag::on_reverse_strand)) {
        flag |= seqan3::sam_flag::on_reverse_strand;
    }
    segment.flag() = flag;
    segment.reference_id() = splitrecord.reference_id();
    segment.reference_position() = refPos;
    segment.cigar_sequence() = cigar;
    segment.sequence() = seq;
    segment.base_qualities() = qual;
//    segment.mapping_quality() = splitrecord.mapping_quality();
    segment.tags() = tags;
    curated.push_back(segment);
}

void SplitReadCalling::sort(const std::string& inputFile, const std::string& outputFile) {
    const char *input = inputFile.c_str();
    const char *output = outputFile.c_str();

    samFile *in = sam_open(input, "rb"); // open input BAM file
    if (!in) {
        std::cerr << helper::getTime() << "Failed to open input BAM file " << inputFile << "\n";
        exit(EXIT_FAILURE);
    }
    samFile *out = sam_open(output, "wb");
    if (!out) {
        std::cerr << helper::getTime() << "Failed to open output BAM file: " << outputFile << "\n";
        sam_close((htsFile *) input);
        exit(EXIT_FAILURE);
    }

    bam_hdr_t *header = sam_hdr_read(in); // load the header
    if (!header) {
        std::cerr << helper::getTime() << "Failed to read header from " << inputFile << "\n";
        sam_close((htsFile *) input);
        exit(EXIT_FAILURE);
    }
    if (sam_hdr_write(out, header) < 0) { // write the header
        std::cerr << helper::getTime() << "Error: Failed to write header to output file " << outputFile << std::endl;
        bam_hdr_destroy(header);
        sam_close(in);
        sam_close(out);
        exit(EXIT_FAILURE);
    }

    // Read BAM records into a vector
    std::vector<bam1_t*> records;
    bam1_t* b = bam_init1();
    while (sam_read1(in, header, b) >= 0) {
        bam1_t* new_record = bam_dup1(b);
        records.push_back(new_record);
    }
    bam_destroy1(b);
    sam_close(in); // close the input file

    // sort the records by QNAME
    std::sort(records.begin(), records.end(), [](const bam1_t* a, const bam1_t* b) {
        return strcmp(bam_get_qname(a), bam_get_qname(b)) < 0;
    });

    // Write sorted BAM records to output file
    for (auto rec : records) {
        if (sam_write1(out, header, rec) < 0) {
            std::cerr << "Error writing BAM record to output file" << std::endl;
            break;
        }
        bam_destroy1(rec);
    }
    // Clean up
    bam_hdr_destroy(header);
    sam_close(out);
}



TracebackResult SplitReadCalling::complementarity(dtp::DNASpan& seq1, dtp::DNASpan& seq2) {
    std::string seq1_str;
    std::string seq2_str;

    if(seq1.size() == 0 || seq2.size() == 0) { return {"","",0,0,0,-1.0,0.0}; }
    for(unsigned z=seq1.size();z-- > 0;) { seq1_str += seq1[z].to_char(); } // rna1 (reverse)
    for(unsigned y=0;y<seq2.size();++y) { seq2_str += seq2[y].to_char(); } // rna2

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

std::pair<double,std::string> SplitReadCalling::hybridization(std::span<seqan3::dna5> &seq1, std::span<seqan3::dna5> &seq2) {
    double hybScore = 1000.0;
    if(seq1.empty() || seq2.empty()) { return {1000,""}; }
    std::string rna1str, rna2str = "";
    for(unsigned i=0;i<seq1.size();++i) { rna1str += seq1[i].to_char(); } // convert to string
    for(unsigned i=0;i<seq2.size();++i) { rna2str += seq2[i].to_char(); }

    // remove non-printable
    rna1str = helper::removeNonPrintable(rna1str); // remove non printable
    rna2str = helper::removeNonPrintable(rna2str); // remove non printable

    // remove non-ATGC
    rna1str = seqIO::removeNonATGC(rna1str); // remove non ATGC
    rna2str = seqIO::removeNonATGC(rna2str); // remove non ATGC

    std::string hyb = "echo '" + rna1str + "&" + rna2str + "' | RNAcofold"; // call

    FILE* pipe = popen(hyb.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error opening pipe" << std::endl;
        return {1000,""};
    }
    char buffer[1024];
    std::string result;
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result += buffer;
    }
    pclose(pipe);
    if(result.empty()) { return {1000,""}; }

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
    hybScore = stod(matches[0]);

    return std::make_pair(hybScore, dotbracket);
}

void SplitReadCalling::addComplementarityToSamRecord(SAMrecord &rec1, SAMrecord &rec2, TracebackResult &res) {
    rec1.tags().get<"XM"_tag>() = res.matches; // number of matches
    rec2.tags().get<"XM"_tag>() = res.matches;
    rec1.tags().get<"XL"_tag>() = res.length; // length of alignment
    rec2.tags().get<"XL"_tag>() = res.length;
    rec1.tags().get<"XC"_tag>() = (float)res.cmpl; // complementarity
    rec2.tags().get<"XC"_tag>() = (float)res.cmpl;
    rec1.tags().get<"XR"_tag>() = (float)res.ratio; // sitelenratio
    rec2.tags().get<"XR"_tag>() = (float)res.ratio;
    rec1.tags().get<"XA"_tag>() = res.a; // alignment
    rec2.tags().get<"XA"_tag>() = res.b;
    rec1.tags().get<"XS"_tag>() = res.score; // score
    rec2.tags().get<"XS"_tag>() = res.score;
}

void SplitReadCalling::addHybEnergyToSamRecord(SAMrecord &rec1, SAMrecord &rec2, std::pair<double, std::string> &hyb) {
    rec1.tags().get<"XE"_tag>() = (float)hyb.first;
    rec2.tags().get<"XE"_tag>() = (float)hyb.first;
    rec1.tags().get<"XD"_tag>() = hyb.second;
    rec2.tags().get<"XD"_tag>() = hyb.second;
}

void SplitReadCalling::writeStats() {
    if(params["stats"].as<std::bitset<1>>() == std::bitset<1>("1")) {
        fs::path outdir = fs::path(params["outdir"].as<std::string>());
        stats->writeStats(outdir, "detect");
    }
}


