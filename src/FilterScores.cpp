#include "FilterScores.hpp"

struct complMatch : public std::true_type {
    int operator()(seqan3::dna5_vector::value_type a, seqan3::dna5_vector::value_type b) {
        if(a == seqan3::complement(b)) {
            return 1;
        } else {
            if(a == 'G'_dna5 && b == 'T'_dna5 || a == 'T'_dna5 && b == 'G'_dna5) {
                return 1;
            } else {
                return -1;
            }
        }
    }
};

FilterScores::FilterScores() {
    // initialize scoring scheme
    this->scheme = seqan3::align_cfg::gap_cost_affine{
            seqan3::align_cfg::open_score{-3},
            seqan3::align_cfg::extension_score{-2}
    };

}
FilterScores::~FilterScores() {}

void FilterScores::computeComplementarity(dtp::DNASpan seq1, dtp::DNASpan seq2) {
}




