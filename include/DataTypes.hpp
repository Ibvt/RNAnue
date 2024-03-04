#ifndef RNANUE_DATATYPES_HPP
#define RNANUE_DATATYPES_HPP

// boost
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

namespace dtp {
    using PathVector = std::vector<fs::path>;
    using SubPathsMap = std::map<std::string, PathVector>;

    using DNAVector = std::vector<seqan3::dna5>; // seqan3 also provides dna5_vector
    using DNASpan = std::span<seqan3::dna5>;
    using QualSpan = std::span<seqan3::phred42>;

    using QualVector = std::vector<seqan3::phred42>;

    // FASTQ
    using FASTQFormat = seqan3::type_list<std::string, DNASpan, QualSpan>;
    using FASTQFields = seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual>;
    using FASTQRecord = seqan3::record<FASTQFormat, FASTQFields>;

    // data types used in preprocessing (state transition table)
    using Left = std::size_t; // left matching block
    using Right = std::pair<std::size_t, std::size_t>; // right matching block

    using StateTransitionTable = std::map<std::pair<int,seqan3::dna5>,std::tuple<int,int,int,int>>;
    using Tables = std::map<std::pair<std::string,DNAVector>, StateTransitionTable>;
    using State = std::tuple<std::string, int, std::pair<int, int>, std::size_t>;

    using Bases = std::map<std::pair<std::string,DNAVector>,DNAVector>;

    using STTEntry = std::tuple<int,int,int,int>;


}



#endif //RNANUE_DATATYPES_HPP
