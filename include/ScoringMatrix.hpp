#ifndef RNANUE_SCORINGMATRIX_HPP
#define RNANUE_SCORINGMATRIX_HPP

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <map>
#include <bitset>
#include <algorithm>

typedef struct {
    size_t width;
    size_t height;
    int *scoring; // scores
    std::bitset<3> *traceback; //
} ScoringMatrix;

/* traceback matrix
 * 0 indicates diagonal
 * 1 indicates left
 * 2 indicates right
 */

ScoringMatrix createScoringMatrix(const char *a,
                                  const char *b,
                                  int matchScore,
                                  int mismatchScore,
                                  int gapCost);

void freeMatrix(ScoringMatrix *mat);
void printMatrix(ScoringMatrix mat);

#endif //RNANUE_SCORINGMATRIX_HPP