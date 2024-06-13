#ifndef RNANUE_DATATYPES_HPP
#define RNANUE_DATATYPES_HPP

// boost
#include <boost/filesystem.hpp>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>


namespace fs = boost::filesystem;

namespace dtp {
    using PathVector = std::vector<fs::path>;
    using SubPathsMap = std::map<std::string, PathVector>;

    using DNAVector = std::vector<seqan3::dna5>; // seqan3 also provides dna5_vector
    using DNASpan = std::span<seqan3::dna5>;
    using QualSpan = std::span<seqan3::phred42>;
    using QualVector = std::vector<seqan3::phred42>;
    using CigarVector = std::vector<seqan3::cigar>;

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

    // FeaturesFields for GFF3/GTF
    struct FeatureFields {
        std::string seqid;
        std::string source;
        std::string type;
        int start;
        int end;
        std::string score;
        char strand;
        char phase;
        std::string attributes;
        FeatureFields() : seqid(""), source(""), type(""), start(-1), end(-1), score(""),
            strand(' '), phase(' '), attributes("") {}
    };

    struct Feature {
        int start;
        int end;
        char strand;
        std::string id;
        std::string name;
        std::string biotype;
        Feature() : start(-1), end(-1) {}
        Feature(int start, int end, char strand, std::string id, std::string name, std::string biotype) :
            start(start), end(end), strand(strand), id(id), name(name), biotype(biotype) {}
    };

    // IBPTree
    using Interval = std::pair<size_t,size_t>;

    // Stats
    struct StatsFields {
        int readsCount;
        int alignedCount;
        int splitsCount;
        int multSplitsCount;
        int nSurvivedCount;
    };
    using StatsMap = std::map<std::string, StatsFields>;
    using SpliceJunctions = std::map<std::string, std::vector<std::pair<size_t,size_t>>>;

    // Analysis
    using PdfMap = std::map<std::pair<std::string, std::string>, double>;
    struct IntPartner {
        std::string partner;
        std::string orientation;
        IntPartner() : partner(""), orientation("") {}
        IntPartner(std::string partner, std::string orientation) : partner(partner), orientation(orientation) {}

        bool operator<(const IntPartner& other) const {
            if(partner < other.partner) {
                return true;
            } else if(partner == other.partner) {
                return orientation < other.orientation;
            } else {
                return false;
            }
        }
    };
    using IntKey = std::vector<IntPartner>;
}

#endif //RNANUE_DATATYPES_HPP
