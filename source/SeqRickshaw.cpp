#include "SeqRickshaw.hpp"

SeqRickshaw::SeqRickshaw(po::variables_map _params) {
    std::cout << "created instance of sequence RickShaw" << std::endl;
    
    // retrieve parametes for the preprocessing
    minlen = _params["minlen"].as<int>();
    phred = _params["quality"].as<int>();

    // which adapters to trim (5'/3'/both)
    modus = _params["modetrm"].as<int>();

    // readtype (SE or PE)
    readtype = _params["readtype"].as<std::string>();

    std::map<std::string,std::string> adpts;

    std::string adpt5File = _params["adpt5"].as<std::string>();
    std::string adpt3File = _params["adpt3"].as<std::string>();

    // create folder for outputs of Rickshaw (e.g., lookup tables)
    fs::path outDirSub{_params["outdir"].as<std::string>()};
    outDirSub /= fs::path(_params["subcall"].as<std::string>());

    try {
        switch(modus) {
            case 0: adpt5Table = calcLookupTable("5\'-adapter", adpt5File);
                    break;
            case 1: 
            {
                adpt3Table = calcLookupTable("3\'-adapter", adpt3File);
                fs::path outDirSub3Adpt = outDirSub / fs::path("table_3adpt.txt");
                std::ofstream lookup3AdptOut(outDirSub3Adpt.string(), std::ios::trunc);
                writeLookupTable(lookup3AdptOut);
                break;
            }
            case 2: adpt5Table = calcLookupTable("5\'-adapter", adpt5File);
                    adpt3Table = calcLookupTable("3\'-adapter", adpt3File);
        }
    }
    catch (seqan3::file_open_error& err) {
        std::cout << err.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }

}

// calculate the smart state transition table for the patterns
std::map<std::pair<std::string,std::string>, LookupTable> SeqRickshaw::calcLookupTable(std::string _type, std::string _adptFile) {
    std::cout << "create lookup table using " << _adptFile << std::endl;

    // 
    Adapters tables;
    seqan3::sequence_file_input adpt{_adptFile}; // read the adapters file
    
    using record_type = typename decltype(adpt)::record_type;
    std::vector<record_type> records{};
//    std::ranges::copy(adpt, std::ranges::back_inserter(records)); -> didn't work when compiling on tesla

    for(auto& rec : adpt) {
        auto readID = seqan3::get<seqan3::field::id>(rec) | seqan3::views::to_char;
        auto readSeq = seqan3::get<seqan3::field::seq>(rec) | seqan3::views::to_char;

        std::pair<std::string, std::string> readAsStrings;
        readAsStrings = std::make_pair(std::string(readID.begin(),readID.end()),std::string(readSeq.begin(),readSeq.end()));

        LookupTable bla = calcShift(readSeq);
        tables.insert(std::make_pair(readAsStrings, bla));
    }
    return tables;
}


/*
 *  
 */
LookupTable SeqRickshaw::calcShift(auto sequence) {
    typedef std::tuple<long,long,long,std::bitset<1>> Entry;

    LookupTable lookup;

    // determine alphabet of sequence
    std::vector<char> alphabet = determineAlphabet(sequence);

        
    // intitial state - reading position to the right (e.g., 000000*)
    // create states set and add intial state
    // std::string state = std::string(sequence.size()-1,'0')+'*';
    std::string state = std::string(sequence.size(),'0');
    States states;

        std::size_t left = std::string::npos; // number of chars in the left block
        std::pair<std::size_t,std::size_t> right = std::make_pair(std::string::npos,std::string::npos);
        size_t readPos = sequence.size()-1;
        
        states.push_back(std::make_tuple(state, left, right,readPos)); //left, right maintain info of previously matched chars 

        //
        int shift = 0;
        int match = 0; // match has been found
        
        std::string nextState;
        int nextStateID = 0;
        int nextReadPosVar = -1;

        int matched = 0; // how many 
        size_t nrMatches = 0;

        int occ;

        std::string suffix;
        int suffixPos;

        char mismatch = ' ';

        std::string pattern;
        std::string subpat;

        size_t matchedCharsCount = 0;

        int goodSuffix = 0;
        int badChar = 0;

        seqan3::debug_stream << sequence << std::endl;   

        States::size_type statesSize = states.size();
        for(States::size_type i=0;i<statesSize;++i) {
            state = std::get<0>(states[i]); // state (e.g., 000X000)
            left = std::get<1>(states[i]); // 
            right = std::get<2>(states[i]);
            readPos = std::get<3>(states[i]);

            std::cout << "-----------------------------------" << std::endl;
            std::cout << sequence << " length: " << sequence.size() << std::endl;
            std::cout << "current state: " << state << std::endl;
            std::cout << "readPos: " << readPos << std::endl;
            std::cout << "left: " << left << std::endl;
            std::cout << "rightStart: " << right.first << std::endl;
            std::cout << "rightEnd: " << right.second << std::endl;
           
            /*
            readPos = calcReadPos(sequence,left,right); // determines the readPos
            if(readPos == std::string::npos) {
                continue;
            }*/

            for(unsigned j=0;j<alphabet.size();++j) { // iterate through any other possibility
                std::cout << "\t---------" << std::endl;
                std::cout << "\talphabet " << alphabet[j] << std::endl;
                std::cout << "\tletter " << sequence[readPos] << " (" << readPos << ")" << std::endl;

                std::size_t nextLeft = left;
                std::pair<std::size_t,std::size_t> nextRight = right;
                nextState = state;
                std::size_t nextReadPos;
                
                if(sequence[readPos] == alphabet[j]) { //  match
                    std::cout << "\tmatch" << std::endl;

                    if(right.first == std::string::npos && right.second == std::string::npos) {
                        nextRight.first = sequence.size()-1;
                        nextRight.second = sequence.size()-1;
                    } else {
                        if(readPos > right.second) {
                            ++nextRight.second;
                        } else {
                            --nextRight.first;
                        }
                    }
                    nextState[readPos] = 'X';
                    std::cout << "\tNextLeft: " << nextLeft << std::endl;
                    std::cout << "\tNextRightStart " << nextRight.first << std::endl;
                    std::cout << "\tNextRightEnd " << nextRight.second << std::endl;
                    std::cout << "\tNextState: " << nextState << std::endl;
                        
                    nextReadPos = calcReadPos(sequence,nextLeft,nextRight); // determines the readPos
                    std::cout << "\tNextReadPos: " << nextReadPos << std::endl;

                    // exclude state with all matches
                    if(std::count(nextState.begin(), nextState.end(), 'X') != sequence.size()) {
                        
                        nextStateID = addState(states, std::make_tuple(nextState, nextLeft, nextRight, nextReadPos), statesSize);
                        std::cout << "\tnextStateID " << nextStateID << '\n';

                        // add to lookup table
                        lookup.insert(std::make_pair(std::make_pair(i, alphabet[j]), std::make_tuple(0,nextStateID,nextReadPos,0)));

                    } else { // a complete match could have been detected
                        // add to lookup table
                        lookup.insert(std::make_pair(std::pair(i, alphabet[j]), std::make_tuple(0,0,0,1)));
                    }


                } else {
                    std::cout << "\tmismatch" << std::endl;

                    // convert range to container
                    auto tmp = sequence | ranges::to<std::vector<char>>;
                    std::string pattern(tmp.begin(),tmp.end());

                    // determine the suffix to check against the pattern
                    if(right.first == std::string::npos && right.second == std::string::npos) {
                        suffix = alphabet[j]; // suffix consists only of mismatch
                    } else {
                        suffix = pattern.substr(right.first, (right.second-right.first)+1);
                        if(readPos < right.first) {
                            suffix = alphabet[j] + suffix;
                        } else {
                            suffix = suffix + alphabet[j];
                        }
                    }

                    std::cout << "\tsuffix " << suffix << std::endl;
                    shift = transition(pattern, suffix, readPos, nextLeft, nextRight);
                    std::cout << "\tLeft " << left << std::endl;
                    
                    std::cout << "\tshift " << shift << std::endl;
                    std::cout << "\tnextLeft " << nextLeft << std::endl;
                    std::cout << "\tnextRightStart " << nextRight.first << std::endl;
                    std::cout << "\tnextRightEnd " << nextRight.second << std::endl;

                    nextState = std::string(sequence.size(),'0');

                    // change stateID
                    if(nextRight.first != std::string::npos && nextRight.second != std::string::npos) {
                        for(unsigned s=nextRight.first;s<=nextRight.second;++s) {
                            nextState[s] = 'X';
                        }
                    }

                    if(nextLeft != std::string::npos) {
                        for(unsigned l=0;l<=nextLeft;++l) {
                            nextState[l] = 'X';
                        }
                    }
                        
                    nextReadPos = calcReadPos(sequence,nextLeft,nextRight); // determines the readPos

                    std::cout << "\tnextReadPos " << nextReadPos << '\n';

                    std::cout << "\tnextState " << nextState << std::endl;
                    nextStateID = addState(states, std::make_tuple(nextState, nextLeft, nextRight, nextReadPos), statesSize);

                    std::cout << "\tnextStateID " << nextStateID << '\n';
                        
                    lookup.insert(std::make_pair(std::make_pair(i, alphabet[j]), std::make_tuple(shift, nextStateID, nextReadPos, 0)));
                }
            }
        }
        std::cout << "list all states" << std::endl;
        States::iterator it;
        for(it = states.begin();it!=states.end();++it) {
            std::cout << std::get<0>(*it) << "\tLeft:" << std::get<1>(*it) << "\tRightStart:" << std::get<2>(*it).first << "\tRightEnd:" << std::get<2>(*it).second;
            std::cout << "\treadPos: " << std::get<3>(*it) << std::endl;
        }

        std::cout << "lookup table" << std::endl;
        for(auto const& x: lookup) {
            std::cout << "(" << x.first.first << "," << x.first.second << ") -> ";
            std::cout << "(" << std::get<0>(x.second) << "," << std::get<1>(x.second) << "," << std::get<2>(x.second) <<  "," << std::get<3>(x.second) << ")" << std::endl;
        }

    return lookup;

}
        
std::size_t SeqRickshaw::calcReadPos(auto& sequence, std::size_t& left, std::pair<std::size_t,std::size_t>& right) {
    std::size_t readPos;

//    std::cout << "readPos: "; 

    // determine readPos
    if(right.first == std::string::npos && right.second == std::string::npos) {
        readPos = sequence.size()-1; // start with first element
 //       std::cout << readPos << " - start with first element (to the right)" << std::endl;
    } else {
        // check if all positions are matches 
        if((right.second) - right.first == sequence.size() - 1 
                || left + (right.second - right.first) == sequence.size() - 2) {
  //          std::cout << "all matched " << std::endl;
            return std::string::npos;
//            continue;
        } else {
            if(right.second + 1 < sequence.size()) { // right end has not been reached
                readPos = right.second + 1; // go to the right
   //             std::cout << readPos << " - go to the right (there is still space) " << std::endl;
            } else { // continue to the left
                if(left == std::string::npos || left < right.first) {
                    readPos = right.first - 1;
    //                std::cout << readPos << " - go to the left " << std::endl;
                }
            } 
        }
    }

    return readPos;
}


        
int SeqRickshaw::transition(std::string pattern, std::string suffix, int readPos, std::size_t& left, std::pair<std::size_t,std::size_t>& right) {
    int shift = 0;


    std::string subsuffix;


    std::size_t localShift = 0;

    // determine left most position of suffix + mismatch
    int leftMostPos = right.first;
    if(right.first == std::string::npos && right.second == std::string::npos) {
        leftMostPos = pattern.size()-1;
    } else {
        if(readPos < right.first) {
            leftMostPos = right.first - 1;
        }
    }

    std::vector<std::size_t> occs;
    findAllOcc(occs, pattern, suffix);

    std::size_t fnd = std::string::npos;

    //
    for(unsigned i=0;i<occs.size();++i) {
        localShift = leftMostPos - occs[i];
        if(left == std::string::npos || localShift - 1 >= left) {
            fnd = occs[i];
            break;
        }
    }

    // determine shift and the blocks (left, right)
    if(fnd != std::string::npos) {
        if(fnd != 0) {
            left = std::string::npos;
            right.first = fnd;
            right.second = fnd + (suffix.size()-1);
            shift = leftMostPos - fnd;
        } else {
            left = suffix.size()-1;   
            right.first = std::string::npos;
            right.second = std::string::npos;
            shift = leftMostPos - fnd;
        }
    } else {
        for(unsigned j=1;j<suffix.size();++j) {
            subsuffix = suffix.substr(j,suffix.size()-j);
            std::cout << "\t\tsubsuffix: " << subsuffix << std::endl;

            if(pattern.substr(0,subsuffix.size()) == subsuffix) {
                fnd = 0;
                std::cout << "\t\tfound " << std::endl;

                left = subsuffix.size()-1;
                right.first = std::string::npos;
                right.second = std::string::npos;
                shift = leftMostPos + j;
            } 
        }
        if(fnd == std::string::npos) {
            left = std::string::npos;
            right.first = std::string::npos;
            right.second = std::string::npos;
            shift = suffix.size();
        }
    }

    return shift;
}



/*
 * return the position within vector
 * */
int SeqRickshaw::addState(States &states, State state, States::size_type &size) {
    // only add state if not in states
    int pos = -1;
    States::iterator it = std::find(states.begin(), states.end(), state);
    if(it != states.end()) { // state is already within states
        pos = std::distance(states.begin(), it);
    } else { // state is not (yet) in states
        states.push_back(state); // push state to vector
        // stateIDs starting from 0
        pos = ++size - 1; // 
    }
    return pos;
}


// determine next read position in state
int SeqRickshaw::nextReadPos(std::string state, int currReadPos) {
    int found = -1;
    // check if readPos is not on the first
    if(currReadPos != state.size()-1) { 
        // search to the right of the current read position
        for(unsigned j=currReadPos+1;j<state.size();++j) {
            if(state[j] == '0') { // found first 
                std::cout << "\tfound in right at position: " << found <<  std::endl;
                found = j;
                return found;
            }
        }
    }

    // continue searching in the left part
    if(currReadPos != 0) {
        for(unsigned k=currReadPos-1;k>=0;--k) {
            if(state[k] == '0') {
                found = k;
                std::cout << "\tfound in left at position: " <<  found << std::endl;
                return found;
            }
        }
    }

    return -1;
}


/*
    auto bla =  seqan3::get<seqan3::field::seq>(_records[0]);
    auto dingens = bla | seqan3::views::to_char;
    std::cout << dingens.size() << std::endl;
//    seqan3::debug_stream << (_records | seqan3::views::get<0> | seqan3::views::to_char) << std::endl; */


/*
 * TODO: use action unique to implement this function more conventiently in c++20
 * */
std::vector<char> SeqRickshaw::determineAlphabet(auto _seq) {
    std::vector<char> alphabet{};
    for(unsigned i=0; i < _seq.size(); ++i) {
        if(std::find(alphabet.begin(), alphabet.end(), _seq[i]) == alphabet.end()) {
            alphabet.push_back(_seq[i]);
        }
    }
    return alphabet;
}

void SeqRickshaw::findAllOcc(std::vector<std::size_t>& fnd, std::string str, std::string substr) {
    std::size_t pos = str.find(substr);
    while(pos != std::string::npos) {
        fnd.push_back(pos);

        pos = str.find(substr, pos + substr.size());
    }
    std::reverse(fnd.begin(),fnd.end());
}



// the boyer moore search algorithm
std::size_t SeqRickshaw::boyermoore(auto& read, LookupTable tab, int m) {
    std::size_t align = 0; 
    int state = 0;
    int readPos = m - 1; // length of pattern to match against
    std::size_t shift = 0;


    std::tuple<int,int,int,int> entry;


    while(align < read.size() - m) {
        auto c = (read | seqan3::views::to_char)[align + readPos];

        if(tab.find(std::make_pair(state,c)) == tab.end()) {
            shift = 1;
            state = 0;
            readPos = m - 1;

        } else {
            entry = tab[std::make_pair(state,c)];
            shift = std::get<0>(entry);
            state = std::get<1>(entry); 
            readPos = std::get<2>(entry);
        
            if(std::get<3>(entry) == 1) {
                return align;
            }
        }
        align = align + shift;
    }
    return std::string::npos;
}


//
std::pair<std::size_t,std::size_t> SeqRickshaw::trimming(auto& fwd) {
    // contains the search positons
    std::size_t fndPos;
    std::vector<size_t> res3Adpt;
    std::vector<size_t> res5Adpt;

    // marks the boundaries
    std::pair<std::size_t,std::size_t> bnds = std::make_pair(0,fwd.size());

    // trim 5' adapters
    if(modus == 0 || modus == 2) {
        for(auto const& x : adpt5Table) {
            fndPos = boyermoore(fwd, x.second, x.first.second.size());
            if(fndPos == 0) { // adapter needs to be detected beginning with the first
                bnds.first = x.first.second.size();
            }
        }
    }

    // trim 3' adapters
    if(modus == 1 || modus == 2) {
        for(auto const& x : adpt3Table) {
            fndPos = boyermoore(fwd, x.second, x.first.second.size());
            if(fndPos < bnds.second) {
                bnds.second = fndPos;
            }
        }
    }
    return bnds;
}



SeqRickshaw::SeqRickshaw() {
}

void SeqRickshaw::start(pt::ptree sample) {
    std::cout << "start" << std::endl;
    
    std::pair<std::size_t,std::size_t> bndsFwd; //
    std::pair<std::size_t,std::size_t> bndsRev; //

    pt::ptree input = sample.get_child("input");
    pt::ptree output = sample.get_child("output");

    std::string forward = input.get<std::string>("forward");
    seqan3::sequence_file_input fwd{forward};
                

    //
    if(readtype == "SE") {
        for(auto & [seq, id, qual] : fwd) {
            bndsFwd = trimming(seq);
        }
    } else { // readtype == "PE"
        std::string reverse = input.get<std::string>("reverse");
        seqan3::sequence_file_input rev{reverse};
   
        // create output files for r1only/r2only - reads that don't survive the trimming
        
        std::string r1only = output.get<std::string>("R1only");
        std::string r2only = output.get<std::string>("R2only");
        seqan3::sequence_file_output r1onlyOut{r1only};
        seqan3::sequence_file_output r2onlyOut{r2only};
        
        for(auto && [rec1,rec2] : seqan3::views::zip(fwd,rev)) {

            bndsFwd = trimming(seqan3::get<seqan3::field::seq>(rec1));
            auto trmReadFwdID = seqan3::get<seqan3::field::id>(rec1);
            auto trmReadFwd = seqan3::get<seqan3::field::seq>(rec1) | seqan3::views::slice(bndsFwd.first,bndsFwd.second);
            auto trmReadFwdQual = seqan3::get<seqan3::field::qual>(rec1) | seqan3::views::slice(bndsFwd.first,bndsFwd.second);
            bool filtFwd = filtering(rec1) && (std::ranges::size(trmReadFwd) >= minlen);

            bndsRev = trimming(seqan3::get<seqan3::field::seq>(rec2));
            auto trmReadRevID = seqan3::get<seqan3::field::id>(rec2);
            auto trmReadRev = seqan3::get<seqan3::field::seq>(rec2) | seqan3::views::slice(bndsRev.first,bndsRev.second);
            auto trmReadRevQual = seqan3::get<seqan3::field::qual>(rec2) | seqan3::views::slice(bndsRev.first,bndsRev.second);
            bool filtRev = filtering(rec2) && (std::ranges::size(trmReadRev) >= minlen);

            if(filtFwd && filtRev) {
                merging(trmReadFwd,trmReadRev);


                // perform merging operation
                /*
                std::cout << trmReadFwdID << std::endl;
                std::cout << (trmReadFwd | seqan3::views::to_char) << std::endl;
                std::cout << trmReadRevID << std::endl;
                std::cout << (trmReadRev | seqan3::views::to_char) << std::endl;
                */
            } else {
                if(filtFwd) {
                    // push to r1only
                    r1onlyOut.emplace_back(trmReadFwd,trmReadFwdID,trmReadFwdQual);
                }
                if(filtRev) {
                    // push to r2only
                    r2onlyOut.emplace_back(trmReadRev,trmReadRevID,trmReadRevQual);
                }
            }

            /*
            std::cout << bndsRev.second<< std::endl;
            std::cout << "original read" <<  (seqan3::get<seqan3::field::seq>(rec2) | seqan3::views::to_char) << std::endl;
            std::cout << "trimmed read" <<  (trmReadRev | seqan3::views::to_char) << std::endl;
            */
        }
    }


//    seqan3::sequence_file_input reads{infwd};
 //   for (auto & [seq, id, qual] : reads )
  //  {
        //processing(seq, adpt3Table);
        /*
        seqan3::debug_stream << "ID: "        << id   << '\n';
        seqan3::debug_stream << "SEQ: "       << seq  << '\n';
        seqan3::debug_stream << "EMPTY QUAL." << qual << '\n'; // qual is empty for FastA files
        */
   // }
}


//
std::string SeqRickshaw::longestCommonSubstr(std::string s1, std::string s2) {
    // Find length of both the strings. 
    int m = s1.length(); 
    int n = s2.length(); 
  
    // Variable to store length of longest 
    // common substring. 
    int result = 0; 
  
    // Variable to store ending point of 
    // longest common substring in X. 
    int end; 
  
    // Matrix to store result of two 
    // consecutive rows at a time. 
    int len[2][n]; 
  
    // Variable to represent which row of 
    // matrix is current row. 
    int currRow = 0; 

    // For a particular value of i and j, 
    // len[currRow][j] stores length of longest 
    // common substring in string X[0..i] and Y[0..j]. 
    for (int i = 0; i <= m; i++) { 
        for (int j = 0; j <= n; j++) { 
            if (i == 0 || j == 0) { 
                len[currRow][j] = 0; 
            } 
            else if (s1[i - 1] == s2[j - 1]) { 
                len[currRow][j] = len[1 - currRow][j - 1] + 1; 
                if (len[currRow][j] > result) { 
                    result = len[currRow][j]; 
                    end = i - 1; 
                } 
            } 
            else { 
                len[currRow][j] = 0; 
            } 
        } 
  
        // Make current row as previous row and 
        // previous row as new current row. 
        currRow = 1 - currRow; 
    } 
  
    // If there is no common substring, print -1. 
    if (result == 0) { 
        return "-1"; 
    } 
  
    // Longest common substring is from index 
    // end - result + 1 to index end in X. 
    return s1.substr(end - result + 1, result); 

}





void SeqRickshaw::merging(auto fwd, auto rev) {
    std::cout << "merging " << '\n';

    auto forward = fwd | seqan3::views::to_char;
    auto reverse = rev | seqan3::views::complement | std::views::reverse | seqan3::views::to_char;


    std::cout << (rev | seqan3::views::to_char) << std::endl;
    std::cout << (reverse | seqan3::views::to_char) << std::endl;



    std::string s1(forward.begin(),forward.end());
    std::string s2(reverse.begin(),reverse.end());


    for(unsigned i=1;i<s1.size()/2;++i) {
//        std::string subs = s1.substr(s1

        std::cout << s1.substr(s1.size()-1-i,i) << std::endl;


    }



    
    
    /*
    std::size_t n;
    for(std::size_t i=0;i<n;++i) {
        if(std::strncmp(lef

    }
    */


}


bool SeqRickshaw::filtering(auto& rec) {
    auto qual = seqan3::get<seqan3::field::qual>(rec) | std::views::transform([] (auto q) { return q.to_phred(); });
    double sum = std::accumulate(qual.begin(), qual.end(), 0);
    double phredScore = sum / std::ranges::size(qual);

    if(phredScore >= phred) {
        return true;
    } else {
        return false;
    }
}



// write lookup table
void SeqRickshaw::writeLookupTable(std::ofstream& ofs) {
    for(auto const& [key, val] : adpt3Table) {
        ofs << "ID: " << key.first << '\n';
        ofs << "Seq: " << key.second << '\n';
        for(auto const& [key2, val2] : val) {
            ofs << '(' << key2.first;
            ofs << ',' << key2.second << ") ->";
            ofs << '\t' << std::get<0>(val2);
            ofs << '\t' << std::get<1>(val2);
            ofs << '\t' << std::get<2>(val2);
            ofs << '\t' << std::get<3>(val2) << '\n';
        }
        ofs << '\n';
    }
}


