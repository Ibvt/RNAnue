#include "ScoringMatrix.hpp"

// determine max element of array
std::pair<int,std::bitset<3>> findMaxElements(double scores[], std::size_t arrSize) {
    std::bitset<3> trace("001"); // remember where 
    int maxElement = scores[0];

    for(int i=1;i<arrSize;++i) {
        if(scores[i] > maxElement && scores[i] > 0) {
            maxElement = scores[i];
            trace.reset();
            trace.set(i);
        }
        if(scores[i] == maxElement && scores[i] > 0) {
            trace.set(i);
        }
    }
    if(maxElement <= 0) {
        maxElement = 0;
        trace.reset();
    }
    return std::make_pair(maxElement, trace);
}


// create the scoring matrix
ScoringMatrix createScoringMatrix(const char *a, const char *b, int matchScore, int mismatchScore, int gapCost) {	
	size_t height = strlen(a)+1;
	size_t width = strlen(b)+1;

	int *scoring = (int*)calloc(width * height, sizeof(int));
    std::bitset<3> *traceback = (std::bitset<3>*) malloc(width * height * sizeof(std::bitset<3>));

	ScoringMatrix matrix{width, height, NULL, NULL}; // create matrix object
    
    for(int i=0;i<matrix.height;++i) {
        for(int j=0;j<matrix.width;++j) {
            *(traceback + i * matrix.width + j) = std::bitset<3>("000");
        }
    }
   
	// determine scoring scheme (complementarity)
	std::map<char, std::map<char,int>> scoringScheme;
	const char *alphabet = "ATGC";
    /*
	for(char c = *alphabet; c; c=*++alphabet) {
		for(char d = *alphabet; d; d=*++alphabet) {
			if(((c == 'A') && (d == 'T')) 
				|| (c == 'T' && d == 'A')
				|| (c == 'G' && (d == 'C'|| d == 'T'))
				|| (c == 'C' || (c == 'T' && d == 'G'))) {
				scoringScheme[c][d] = matchScore;
			} else {
				scoringScheme[c][d] = mismatchScore;
			}
		} 
	}*/

    scoringScheme['A']['T'] = matchScore;
    scoringScheme['T']['A'] = matchScore;
    scoringScheme['G']['C'] = matchScore;
    scoringScheme['C']['G'] = matchScore;
    scoringScheme['G']['T'] = matchScore;
    scoringScheme['T']['G'] = matchScore;

    scoringScheme['T']['T'] = mismatchScore;
    scoringScheme['A']['A'] = mismatchScore;
    scoringScheme['C']['C'] = mismatchScore;
    scoringScheme['G']['G'] = mismatchScore;

    scoringScheme['A']['G'] = mismatchScore;
    scoringScheme['A']['C'] = mismatchScore;
    scoringScheme['G']['A'] = mismatchScore;
    scoringScheme['C']['A'] = mismatchScore;

    scoringScheme['T']['C'] = mismatchScore;
    scoringScheme['C']['T'] = mismatchScore;

	if(scoring != NULL && traceback != NULL) { // memory allocation was successful
		for(unsigned i=1; i<height; ++i) {
			for(unsigned j=1; j<width; ++j) {
                double *scores;
                scores = (double*)calloc(3,sizeof(double));
                if(scores) { // memory allocation for scores successful
                    scores[0] = *(scoring + ((i-1) * width) + (j-1)) + scoringScheme[a[i-1]][b[j-1]];
                    scores[1] = *(scoring + ((i-1) * width) + j) - gapCost;
                    scores[2] = *(scoring + (i * width) + (j-1)) - gapCost;

                    /*
                    std::cout << "scores / compare: " << a[i-1] << " and " << b[j-1] << std::endl;
                    std::cout << "similarity: " << scoringScheme[a[i-1]][b[j-1]] << std::endl;
                    std::cout << scores[0] << std::endl;
                    std::cout << scores[1] << std::endl;
                    std::cout << scores[2] << std::endl;
                    std::cout << sc.second << std::endl;
                    */

                    std::pair<int,std::bitset<3>> sc = findMaxElements(scores,3);
                    *(scoring + i * width + j) = sc.first;
                    *(traceback + i * width + j) = sc.second;
                }
			}
		}
        matrix.scoring = scoring;
        matrix.traceback = traceback;
	} else {
		std::cout << "memory allocation for scoring matrix failed" << "\n";
	}

	return matrix;
}

//
void freeMatrix(ScoringMatrix *matrix) {
    free(matrix->scoring);
    free(matrix->traceback);
    matrix->scoring = NULL;
    matrix->traceback = NULL;
}

//
void printMatrix(ScoringMatrix matrix) {
    std::cout << "Scoring Matrix\n";
    for(int i=0;i<matrix.height;++i) {
        for(int j=0;j<matrix.width;++j) {
            printf("%d ", *(matrix.scoring + i * matrix.width + j));
        }
        printf("\n");
    }
    
    std::cout << "Traceback Matrix\n";
    for(int i=0;i<matrix.height;++i) {
        for(int j=0;j<matrix.width;++j) {
            std::cout << *(matrix.traceback + i * matrix.width + j) << " ";
        }
        std::cout << '\n';
    }
}
