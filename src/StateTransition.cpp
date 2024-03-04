#include "StateTransition.hpp"

StateTransition::StateTransition(std::string type, std::string adptFile) {
    std::cout << helper::getTime() << "Create " << type << " lookup table using " << adptFile << "\n";
    tables = preprocessing(type, adptFile);
}

StateTransition::StateTransition() {}
StateTransition::~StateTransition() {}

dtp::Tables StateTransition::getTables() {
    return tables;
}

dtp::Bases StateTransition::getBases() {
    return bases;
}

dtp::Tables StateTransition::preprocessing(std::string type, std::string adptFile) {
    std::cout << helper::getTime() << "Create lookup table using " << adptFile << "\n";
    dtp::Tables adptTables;

    seqan3::sequence_file_input adptSeqFile{adptFile};
    for(auto& rec : adptSeqFile) {
        dtp::StateTransitionTable tab = calcStateTransitionTable(rec.sequence());
        adptTables.insert(std::make_pair(std::make_pair(rec.id(), rec.sequence()), tab));

        // stores the bases for each sequence (read)
        bases.insert(std::make_pair(std::make_pair(rec.id(), rec.sequence()), seqIO::determineAlphabet(rec.sequence())));
    }
    return adptTables;
}

dtp::StateTransitionTable StateTransition::calcStateTransitionTable(dtp::DNAVector seq) {
    dtp::StateTransitionTable lookup;
    // determine alphabet of sequence
    std::vector<seqan3::dna5> alphabet = seqIO::determineAlphabet(seq);

    // define left (single value) and right block (pair) - left
    std::size_t left = std::string::npos; // left block (set to max)
    std::pair<std::size_t,std::size_t> right = std::make_pair(std::string::npos,std::string::npos);
    size_t readPos = seq.size()-1;

    // create initial state (initialized with all 0's)
    std::string state = std::string(seq.size(), '0'); // intial state (only 0s)
    std::vector<dtp::State> states; // vector to store all the states
    states.push_back(std::make_tuple(state, left, right, readPos));

    int shift = 0;
    int match = 0; // if match has been found (=1)

    std::string nextState;
    int nextStateID = 0;
    int nextReadPosVar = -1;

    int matched = 0;
    size_t nrMatches = 0;

    int occ;

    std::vector<seqan3::dna5> suffix;
    std::span<seqan3::dna5> suffixSpan;

    int suffixPos;
    char mismatch = ' ';

    std::string pattern;
    std::string subpat;

    size_t matchedCharsCount = 0;

    int goodSuffix = 0;
    int badChar = 0;

    std::vector<dtp::State>::size_type statesSize = states.size();
    for(std::vector<dtp::State>::size_type i=0;i<statesSize;++i) {
        // get current state (state, left, right, readPos)
        state = std::get<0>(states[i]);
        left = std::get<1>(states[i]);
        right = std::get<2>(states[i]);
        readPos = std::get<3>(states[i]);

        for(unsigned j=0;j<alphabet.size();++j) { // iterate through any other possibility
            std::size_t nextLeft = left;
            std::pair<std::size_t,std::size_t> nextRight = right;
            nextState = state;
            std::size_t nextReadPos;
            if(seq[readPos] == alphabet[j]) { // match
                if(right.first == std::string::npos && right.second == std::string::npos) { // right block is empty
                    // set right block to last element
                    nextRight.first = seq.size()-1;
                    nextRight.second = seq.size()-1;
                } else {
                    if (readPos > right.second) {
                        ++nextRight.second;
                    } else {
                        --nextRight.first;
                    }
                }
                nextState[readPos] = 'X';
                nextReadPos = calcReadPos(seq, nextLeft, nextRight);

                if(std::count(nextState.begin(), nextState.end(), 'X') != seq.size()) {
                    nextStateID = addState(states, std::make_tuple(nextState, nextLeft, nextRight, nextReadPos), statesSize);
                    lookup.insert(std::make_pair(std::make_pair(i, alphabet[j]), std::make_tuple(0,nextStateID,nextReadPos,0)));
                } else {
                    lookup.insert(std::make_pair(std::pair(i, alphabet[j]), std::make_tuple(0,0,0,1)));
                }

            } else { // mismatch
                // determine the suffix to check against the pattern
                if(right.first == std::string::npos && right.second == std::string::npos) {
                    suffix = std::vector<seqan3::dna5>{alphabet[j]};
                } else {
                    suffix = seqIO::spanToVector(seqan3::views::slice(seq,
                                                                       right.first,
                                                                       right.first+(right.second-right.first)+1));
                    if(readPos < right.first) {
                        suffix.insert(suffix.begin(), alphabet[j]);
                    } else {
                        suffix.push_back(alphabet[j]);
                    }
                }

                shift = transition(seq, suffix, readPos, nextLeft, nextRight);
                nextState = std::string(seq.size(), '0');

                // change state
                if(nextRight.first != std::string::npos && nextRight.second != std::string::npos) {
                    for(unsigned k=nextRight.first;k<=nextRight.second;++k) {
                        nextState[k] = 'X';
                    }
                }

                if(nextLeft != std::string::npos) {
                    for(unsigned l=0;l<=nextLeft;++l) {
                        nextState[l] = 'X';
                    }
                }

                nextReadPos = calcReadPos(seq, nextLeft, nextRight);
                nextStateID = addState(states, std::make_tuple(nextState, nextLeft, nextRight, nextReadPos), statesSize);
                lookup.insert(std::make_pair(std::make_pair(i, alphabet[j]), std::make_tuple(shift,nextStateID,nextReadPos,0)));
            }
        }
    }
    return lookup;
}

// calculate the new read position - always starts to the left of the right block
std::size_t StateTransition::calcReadPos(dtp::DNAVector& seq, dtp::Left& left, dtp::Right& right) {
    std::size_t readPos;

    if(right.first == std::string::npos && right.second == std::string::npos) {
        readPos = seq.size()-1; // start with first element
    } else {
        // check if all positions are matches
        if((right.second) - right.first == seq.size() - 1
           || left + (right.second - right.first) == seq.size() - 2) {
            return std::string::npos;
        } else {
            if(right.second + 1 < seq.size()) { // right end has not been reached
                readPos = right.second + 1; // go to the right
            } else { // continue to the left
                if(left == std::string::npos || left < right.first) {
                    readPos = right.first - 1;
                }
            }
        }
    }
    return readPos;
}

int StateTransition::addState(std::vector<dtp::State>& states, dtp::State state, std::vector<dtp::State>::size_type& size) {
    int pos = -1;
    std::vector<dtp::State>::iterator it = std::find(states.begin(), states.end(), state);
    if(it != states.end()) { // state is already within states
        pos = std::distance(states.begin(), it);
    } else { // state is not (yet) in states
        states.push_back(state); // push state to vector
        // stateIDs starting from 0
        pos = ++size - 1; //
    }
    return pos;
}

int StateTransition::transition(dtp::DNAVector seq,  dtp::DNAVector suffix, int readPos, dtp::Left& left, dtp::Right& right) {

    int shift = 0;
    std::size_t localShift = 0;
    std::vector<seqan3::dna5> subSuffix;
    std::vector<seqan3::dna5> prefix;

    //seqan3::debug_stream << "seq: " << seq << "\n";
    //seqan3::debug_stream << "suffix: " << suffix << "\n";

    int leftMostPos = right.first;
    if(right.first == std::string::npos && right.second == std::string::npos) {
        leftMostPos = seq.size()-1;
    } else {
        if(readPos < right.first) {
            leftMostPos = right.first - 1;
        }
    }


    std::vector<std::size_t> matches = searchPattern(seq, suffix);
    std::size_t fnd = std::string::npos;

    for(unsigned i=0;i<matches.size();++i) {
        localShift = leftMostPos-matches[i];
        if(left == std::string::npos || localShift-1 >= left) {
            fnd = matches[i];
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
            subSuffix = seqIO::spanToVector(seqan3::views::slice(seq, j, suffix.size()));
            // compare subset of suffix with beginning of sequence (mismatch falls off to the left)
            prefix = seqIO::spanToVector(seqan3::views::slice(suffix, 0, subSuffix.size()));
            if(prefix == subSuffix) {
                fnd = 0;
                left = subSuffix.size()-1;
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
// search for suffix in sequence (returns all positions)
std::vector<std::size_t> StateTransition::searchPattern(dtp::DNAVector seq, dtp::DNAVector pat) {
    std::vector<std::size_t> occ; // vector to store the positions of the occurrences
    auto it = seq.begin();
    while((it = std::search(it, seq.end(), pat.begin(), pat.end())) != seq.end()) {
        std::size_t pos = std::distance(seq.begin(), it);
        occ.emplace_back(pos);
        ++it;
    }
    std::reverse(occ.begin(), occ.end());
    return occ;
}

void StateTransition::writeStateTransitionTable(std::ofstream& ofs) {
    std::cout << helper::getTime() << "Write the State Transition Table to file" << std::endl;
    for(auto& [key, val] : tables) {
        ofs << "ID: " << key.first << "\n";
        ofs << "Sequence: ";
        seqIO::printDNAVector(key.second, ofs);
        ofs << "(State,Base)\tShift\tNextState\tNextReadPos\tMatch\n";
        for(auto const& [key2, val2] : val) {
            ofs << "(" << key2.first;
            ofs << "," << seqan3::to_char(key2.second) << ") ->";
            ofs << '\t' << std::get<0>(val2);
            ofs << '\t' << std::get<1>(val2);
            ofs << '\t' << std::get<2>(val2);
            ofs << '\t' << std::get<3>(val2) << '\n';
        }
        ofs << "\n";
    }
}
