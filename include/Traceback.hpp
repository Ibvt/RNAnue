#ifndef TRACEBACK_H
#define TRACEBACK_H

#include "ScoringMatrix.hpp"
#include <vector>

typedef struct {
    std::string a;
    std::string b;
    int length;
    int matches;
    int score;
    double cmpl; // complementarity
    double ratio; // sitelenratio
} TracebackResult;

std::vector<TracebackResult> traceback(ScoringMatrix matrix, const char *a, const char *b);
void free_traceback(TracebackResult trc);

#endif // TRACEBACK_H
