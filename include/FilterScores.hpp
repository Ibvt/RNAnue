#ifndef RNANUE_FILTERSCORES_HPP
#define RNANUE_FILTERSCORES_HPP

#include "DataTypes.hpp"

class Complementarity {
    public:
        Complementarity();
        ~Complementarity();

        void compute(dtp::DNAVector seq1, dtp::DNAVector seq2);

    private:
        double score;
        int alignmentLength;
        double siteLengthRatio;
        int matches;
};

class Hybridization {
    public:
        Hybridization();
        ~Hybridization();

        void compute(dtp::DNAVector seq1, dtp::DNAVector seq2);

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

};

#endif //RNANUE_FILTERSCORES_HPP
