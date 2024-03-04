//
// Created by Richard Albin Schaefer on 1/30/24.
//

#ifndef RNANUE_STATETRANSITION_HPP
#define RNANUE_STATETRANSITION_HPP

#include <iostream>
#include <vector>

// seqan3
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "Utility.hpp"
#include "DataTypes.hpp"

class StateTransition {
    public:
        StateTransition(std::string type, std::string adptFile);
        StateTransition();
        ~StateTransition();

        // getter
        dtp::Tables getTables();
        dtp::Bases getBases();

        // preprocess the search patterns
        dtp::Tables preprocessing(std::string type, std::string adptFile);
        dtp::StateTransitionTable calcStateTransitionTable(dtp::DNAVector seq);
        std::size_t calcReadPos(dtp::DNAVector& seq, dtp::Left& left, dtp::Right& right);
        int addState(std::vector<dtp::State>& states, dtp::State state, std::vector<dtp::State>::size_type& size);
        int transition(dtp::DNAVector seq,  dtp::DNAVector suffix, int readPos, dtp::Left& left, dtp::Right& right);
        std::vector<std::size_t> searchPattern(dtp::DNAVector sequence, dtp::DNAVector pattern);

        void writeStateTransitionTable(std::ofstream& ofs);

    private:
        dtp::Tables tables;
        dtp::Bases bases; // stores the bases for a read (ID, Seq)



};

#endif //RNANUE_STATETRANSITION_HPP
