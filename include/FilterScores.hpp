#ifndef RNANUE_FILTERSCORES_HPP
#define RNANUE_FILTERSCORES_HPP

#include "DataTypes.hpp"

// Seqan3
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/pairwise_combine.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3::literals;

class Complementarity {
    public:
        Complementarity();
        ~Complementarity();

        void compute(dtp::DNAVector seq1, dtp::DNAVector seq2);
};

class Hybridization {
    public:
        Hybridization();
        ~Hybridization();

        void compute(dtp::DNASpan seq1, dtp::DNASpan seq2);

    private:
        double score;
        int alignmentLength;
        double siteLengthRatio;
        int matches;
};


class FilterScores {
    public:
        FilterScores();
        ~FilterScores();

        void computeComplementarity(dtp::DNASpan seq1, dtp::DNASpan seq2);

        seqan3::align_cfg::gap_cost_affine scheme;
};

#endif //RNANUE_FILTERSCORES_HPP
