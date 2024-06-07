#include "Preproc.hpp"

// constructor to process single-end reads
Preproc::Preproc(po::variables_map& params, StateTransition& adpt5Tables, StateTransition& adpt3Tables)
    : params(params), adpt5Tables(adpt5Tables), adpt3Tables(adpt3Tables) {

    // get parameters
    wsize = params["wsize"].as<int>();
    quality = params["quality"].as<int>();
    minlen = params["minlen"].as<int>();
    minovlps = params["minovlps"].as<int>();
}

Preproc::~Preproc() {}

// handling for single-end reads
void Preproc::processing(fs::path& in, fs::path& out) {
    std::cout << helper::getTime() << "Trimming on sample: " << in << "\n";

    std::mutex inMutex, outMutex; //
    seqan3::sequence_file_input fin{in.string()};
    std::ofstream fout(out.string());
    //seqan3::sequence_file_output fout{out.string()}; // causes segfault/assertion error

    std::pair<std::size_t, std::size_t> bnds;
    int readcount = 0;

    auto reads = fin | seqan3::views::async_input_buffer(100);
    auto worker = [&]()
    {
        for (auto & record : reads)
        {
            inMutex.lock();
            std::string seqid = record.id();
            dtp::DNAVector seq = record.sequence();
            dtp::QualVector qual = record.base_qualities();
            if(readcount != 0 && readcount % 100000 == 0) {
                std::cout << helper::getTime() << readcount << " reads processed \n";
            }
            readcount++;
            inMutex.unlock();

            bnds = trimReads(seq); //
            // perform window trimming on the data (if activated)
            if(params["wtrim"].as<std::bitset<1>>() == std::bitset<1>("1")) {
                bnds.second = nibble(qual, bnds);
            }
            dtp::DNASpan trimmedSeq = seq | seqan3::views::slice(bnds.first, bnds.second);
            dtp::QualSpan trimmedQual = qual | seqan3::views::slice(bnds.first, bnds.second);
            dtp::FASTQRecord trimmedRecord{seqid, trimmedSeq, trimmedQual};


            /*
            std::cout << "bnds: " << bnds.first << " " << bnds.second << std::endl;
            std::cout << "Size of orig sequence: " << std::ranges::size(seq) << std::endl;
            std::cout << "Size of orig qual: " << std::ranges::size(seq) << std::endl;
            std::cout << "Size of sequence: " << std::ranges::size(trimmedSeq) << std::endl;
            std::cout << "Size of qual: " << std::ranges::size(trimmedQual) << std::endl;
            assert(std::ranges::size(trimmedSeq) == std::ranges::size(trimmedQual));
             */

            if(calcAvgPhred(trimmedQual) >= quality && trimmedSeq.size() >= minlen) {
                //outMutex.lock();
                std::lock_guard<std::mutex> lock(outMutex);
                fout << fastqToString(trimmedRecord) << "\n";
                //fout.push_back(trimmedRecord); // causes segfault/assertion error
                //outMutex.unlock();
            }
        }
    };

    std::vector<std::future<void>> threadObjects;
    for(unsigned i=0;i<params["threads"].as<int>();++i) {
        threadObjects.push_back(std::async(std::launch::async, worker));
    }

    // close files
    fout.close();
}

// handling for paired-end reads
void Preproc::processing(fs::path& inFwd, fs::path& outFwd, fs::path& inRev, fs::path& outR1only,
                         fs::path& outR2only, fs::path& outR1unmerged, fs::path& outR2unmerged) {

    seqan3::sequence_file_input inFwdSeq{inFwd.string()};
    seqan3::sequence_file_input inRevSeq{inRev.string()};

    std::cout << helper::getTime() << "Preprocessing paired-end reads\n";
    std::cout << helper::getTime() << "forward: " << inFwd << "\n";
    std::cout << helper::getTime() << "reverse: " << inRev << "\n";

    std::ofstream outFwdSeq(outFwd.string());
    std::ofstream outR1onlySeq(outR1only.string());
    std::ofstream outR2onlySeq(outR2only.string());
    std::ofstream outR1unmrgSeq(outR1unmerged.string());
    std::ofstream outR2unmrgSeq(outR2unmerged.string());

    /* causes issues on multicore / assertions error
    seqan3::sequence_file_output outFwdSeq{outFwd.string()};
    seqan3::sequence_file_output outR1onlySeq{outR1only.string()};
    seqan3::sequence_file_output outR2onlySeq{outR2only.string()};
    seqan3::sequence_file_output outR1unmrgSeq{outR1unmerged.string()};
    seqan3::sequence_file_output outR2unmrgSeq{outR2unmerged.string()};*/

    std::mutex inMutex, outMutex;
    std::pair<std::size_t, std::size_t> bndsFwd, bndsRev;

    int readcount = 0;
    auto inFwdSeqBuf = inFwdSeq | seqan3::views::async_input_buffer(100);
    auto inRevSeqBuf = inRevSeq | seqan3::views::async_input_buffer(100);
    auto worker = [&]() {
        for (auto && [fwd, rev]: seqan3::views::zip(inFwdSeqBuf, inRevSeqBuf)) {
            std::string seqidFwd, seqidRev;
            dtp::DNAVector seqFwd, seqRev;
            dtp::QualVector qualFwd, qualRev;
            {
                std::lock_guard lock(inMutex);
                if(readcount != 0 && readcount % 100000 == 0) {
                    std::cout << helper::getTime() << "Processed reads: " << (readcount++) << "\n";
                }
                seqidFwd = fwd.id();
                seqidRev = rev.id();
                seqFwd = fwd.sequence();
                seqRev = rev.sequence();
                qualFwd = fwd.base_qualities();
                qualRev = rev.base_qualities();
            }

            // trimming
            bndsFwd = trimReads(seqFwd);
            bndsRev = trimReads(seqRev);
            if (params["wtrim"].as<std::bitset<1>>() == std::bitset<1>("1")) {
                bndsFwd.second = nibble(qualFwd, bndsFwd);
                bndsFwd.second = nibble(qualRev, bndsRev);
            }

            // retrieve trimmed sequences
            auto trmSeqFwd = seqFwd
                             | seqan3::views::slice(bndsFwd.first, bndsFwd.second);
            auto trmQualFwd = qualFwd
                              | seqan3::views::slice(bndsFwd.first, bndsFwd.second);
            auto trmSeqRev = seqRev
                             | seqan3::views::slice(bndsRev.first, bndsRev.second);
            auto trmQualRev = qualRev
                              | seqan3::views::slice(bndsRev.first, bndsRev.second);

            // perform filtering (quality and length)
            bool fwdFilter = filterReads(trmQualFwd);
            bool revFilter = filterReads(trmQualFwd);

            if(fwdFilter && revFilter) {
                std::pair<dtp::DNAVector, dtp::QualVector> merged = mergeReads(trmSeqFwd, trmQualFwd,
                                                                               trmSeqRev, trmQualRev);

                if (merged.first.size() != 0) {
                    {
                        std::lock_guard<std::mutex> lock(outMutex);
                        dtp::FASTQRecord outMerged{seqidFwd, merged.first, merged.second};
                        outFwdSeq << fastqToString(outMerged) << "\n";
                    }
                } else { // reads could not have been merged
                    {
                        std::lock_guard<std::mutex> lock(outMutex);
                        dtp::FASTQRecord outR1unmrgObj{seqidFwd, trmSeqFwd, trmQualFwd};
                        dtp::FASTQRecord outR2unmrgObj{seqidRev, trmSeqRev, trmQualRev};
                        outR1unmrgSeq << fastqToString(outR1unmrgObj) << "\n";
                        outR2unmrgSeq << fastqToString(outR2unmrgObj) << "\n";
                    }
                }
            } else {
                if(fwdFilter && !revFilter) {
                    {
                        std::lock_guard<std::mutex> lock(outMutex);
                        dtp::FASTQRecord outR1onlyObj{seqidFwd, trmSeqFwd, trmQualFwd};
                        outR1onlySeq << fastqToString(outR1onlyObj) << "\n";
                    }
                } else {
                    if (!fwdFilter && revFilter) { // R2 only
                        // reverse complement the sequence and quality
                        auto revSeqRC = trmSeqFwd | std::views::reverse | seqan3::views::complement;
                        auto revQualR = trmQualFwd | std::views::reverse;
                        dtp::DNAVector revSeqRCVec{revSeqRC.begin(), revSeqRC.end()};
                        dtp::QualVector revQualRVec{revQualR.begin(), revQualR.end()};
                        dtp::FASTQRecord outR2onlyObj{seqidFwd, revSeqRCVec, revQualRVec};
                        {
                            std::lock_guard<std::mutex> lock(outMutex);
                            outR2onlySeq << fastqToString(outR2onlyObj) << "\n";
                        }
                    }
                }
            }
        }
    };

    std::vector<std::future<void>> threadObjects;
    for (unsigned i = 0; i < params["threads"].as<int>(); ++i) {
        threadObjects.push_back(std::async(std::launch::async, worker));
    }

    // close
    outFwdSeq.close();
    outR1onlySeq.close();
    outR2onlySeq.close();
    outR1unmrgSeq.close();
    outR2unmrgSeq.close();
}


//
std::pair<size_t, size_t> Preproc::trimReads(dtp::DNAVector& seq) {
    int modetrm = params["modetrm"].as<int>();
    size_t fndPos = std::string::npos;
    // mark the boundaries (left and right) of the sequence (initially the whole sequence)
    std::pair<size_t, size_t> boundaries = std::make_pair(0, seq.size()-1);

    // trim 5'-adapter
    if(modetrm == 0 || modetrm == 2) {
        for(auto& tab : adpt5Tables.getTables()) {
            dtp::DNAVector bases = adpt5Tables.getBases().find(std::make_pair(tab.first.first, tab.first.second))->second;
            fndPos = boyermoore(seq, tab.second, bases, tab.first.second);
            if (fndPos != std::string::npos) { boundaries.first = fndPos + bases.size()-1; }
        }
    }

    // trim 3'-adapter
    dtp::DNAVector bases;
    if(modetrm == 1 || modetrm == 2) {
        for(auto& tab : adpt3Tables.getTables()) {
            bases = adpt3Tables.getBases().find(std::make_pair(tab.first.first, tab.first.second))->second;
            fndPos = boyermoore(seq, tab.second, bases, tab.first.second);
            if (fndPos != std::string::npos) { boundaries.second = fndPos; }
        }
    }
    return boundaries;
}

std::size_t Preproc::boyermoore(dtp::DNAVector& read, dtp::StateTransitionTable& table,
                                dtp::DNAVector bases, const dtp::DNAVector& pat) {
    int patlen = pat.size();
    int mismatch = patlen * params["mmrate"].as<double>();

    // storing transition (shift, nextState and readPos)
    int align = 0, state = 0;
    int readPos = patlen-1;
    int shift = 0;
    int match = 0;

    dtp::STTEntry transition;
    std::mutex mutex;

    int mmCounter = 0;
    while(align < static_cast<int>(read.size()) - patlen) {
        seqan3::dna5 c = read[align + readPos];
        // make sure that the character can be found in the table
        if(table.find(std::make_pair(state, c)) == table.end()) {
            // c is not in the table - reset boyermoore
            state = 0;
            shift = 1;
            readPos = patlen-1;
            match = 0;
        } else { // c is in the table
            shift = std::get<0>(table[std::make_pair(state, c)]);
            state = std::get<1>(table[std::make_pair(state, c)]);
            readPos = std::get<2>(table[std::make_pair(state, c)]);
            match = std::get<3>(table[std::make_pair(state, c)]);

            // allow mismatches
            if(shift != 0 && mmCounter < mismatch-1) {
                // remove current base (only check for the remainder)
                remove(bases.begin(), bases.end(), c);
                for(unsigned i=0; i<bases.size();++i) {
                    if(std::get<0>(table[std::make_pair(state, bases[i])]) == 0) {
                        shift = 0;
                        state = std::get<1>(table[std::make_pair(state, c)]);
                        // if readPos is 1 (second position), then the pattern is found
                        if(readPos == 1) {
                            match = 1;
                            readPos = 0;
                        } else {
                            readPos = std::get<2>(table[std::make_pair(state, c)]);
                            match = std::get<3>(table[std::make_pair(state, c)]);
                        }
                        break;
                    }
                }
                mmCounter++;
            }
            if(match == 1) {
                // found the pattern
                return align;
            }
        }
        align += shift;
    }
    return std::string::npos;
}

// window trimming method
std::size_t Preproc::nibble(dtp::QualVector& qual, std::pair<std::size_t, std::size_t>& bnds) {
    std::size_t threePrimeEnd = bnds.second;
    while(threePrimeEnd - bnds.first >= wsize) {
        dtp::QualSpan windowScores = qual |
                seqan3::views::slice(threePrimeEnd - wsize, threePrimeEnd);
        double windowPhredAvg = calcAvgPhred(windowScores);

        if(windowPhredAvg >= quality) {
            break;
        }
        threePrimeEnd -= wsize;
    }
    return threePrimeEnd;
}


//
std::pair<dtp::DNAVector, dtp::QualVector> Preproc::mergeReads(dtp::DNASpan& fwdSeq, dtp::QualSpan& fwdQual,
                                                               dtp::DNASpan& revSeq, dtp::QualSpan& revQual) {

    // reverse complement the reverse read
    auto revSeqRC = revSeq | std::views::reverse | seqan3::views::complement;
    auto revQualRC = revQual | std::views::reverse;

    /* convert to String - this seems counterintuitive, but working of views on DNAVectors cannot be concatenated
     * and would have to be converted to a DNAVector afterwards */

    // sequence
    auto fwdSeqTrf = fwdSeq | seqan3::views::to_char;
    auto revSeqRCTrf = revSeqRC | seqan3::views::to_char;
    std::string fwdSeqStr{fwdSeqTrf.begin(), fwdSeqTrf.end()};
    std::string revSeqRCStr{revSeqRCTrf.begin(), revSeqRCTrf.end()};

    // quality scores
    auto fwdQualTrf = fwdQual | seqan3::views::to_char;
    auto revQualTrf = revQualRC | seqan3::views::to_char;
    std::string fwdQualStr{fwdQualTrf.begin(), fwdQualTrf.end()};
    std::string revQualRCStr{revQualTrf.begin(), revQualTrf.end()};

    std::string lcs = longestCommonSubseq(fwdSeqStr, revSeqRCStr);

    if(lcs.size() < 1 || lcs.size() < minovlps) {
        return std::make_pair(dtp::DNAVector{}, dtp::QualVector{});
    } else {
        // search for location of lcs in string
        std::size_t fwdLcsPos = fwdSeqStr.find(lcs);
        std::size_t revLcsPos = revSeqRCStr.find(lcs);

        std::string mrgSeqStr = fwdSeqStr.substr(0, fwdLcsPos+lcs.size()) + revSeqRCStr.substr(0,revLcsPos);
        std::string mrgQualStr = fwdQualStr.substr(0, fwdLcsPos+lcs.size()) + revQualRCStr.substr(0,revLcsPos);

        auto mrgSeq = mrgSeqStr | seqan3::views::char_to<seqan3::dna5>;
        dtp::DNAVector mrgSeqVec{mrgSeq.begin(), mrgSeq.end()};
        auto mrgQual = mrgQualStr | seqan3::views::char_to<seqan3::phred42>;
        dtp::QualVector mrgQualVec{mrgQual.begin(), mrgQual.end()};

        return std::make_pair(mrgSeqVec, mrgQualVec);

    }
}

// calculate the longest common subsequence (LCS) to determine the overlaps
// Note: rev is a transformed object (reverse complement)
std::string Preproc::longestCommonSubseq(std::string& fwd, std::string& rev) {
    // determine lengths
    size_t m = fwd.size();
    size_t n = rev.size();

    int dp[2][n+1]; // dynamic programming table
    int curr=0, res=0, end=0;

    for(unsigned int i=0;i<=m;++i) {
        for(unsigned int j=0;j<=n;++j) {
            if(i==0 || j==0) {
                dp[curr][j]=0;
            } else {
                if(fwd[i-1] == rev[j-1]) {
                    dp[curr][j] = dp[1-curr][j-1]+1;
                    if(res<dp[curr][j]) {
                        res=dp[curr][j];
                        end=i-1;
                    }
                } else {
                    dp[curr][j]=0;
                }
            }
        }
        curr=1-curr;
    }
    if(res==0) {
        return std::string{};
    } else {
        std::string ans;
        return std::string(fwd.substr(end-res+1, res));
    }
}

// filtering reads based on quality and length
bool Preproc::filterReads(auto& qual) {
    // determine the length
    size_t len = std::ranges::size(qual);
    if(calcAvgPhred(qual) >= quality && len >= minlen) {
        return true;
    } else {
        return false;
    }
}

// calculate the average phred score
double Preproc::calcAvgPhred(auto& qual) {
    auto phredScores = qual | std::views::transform([] (auto qual) {
        return seqan3::to_phred(qual);
    });
    auto phredSum = std::accumulate(phredScores.begin(), phredScores.end(), 0);\
    double phredAvg = phredSum / std::ranges::size(phredScores);
    return phredAvg;
}

std::string Preproc::fastqToString(dtp::FASTQRecord& rec) {
    std::stringstream ss;
    ss << "@" << std::get<0>(rec) << "\n";
    for(auto& c : std::get<1>(rec)) { ss << c.to_char(); }
    ss << "\n+\n";
    for(auto& q : std::get<2>(rec)) { ss << q.to_char(); }
    return ss.str();
}
