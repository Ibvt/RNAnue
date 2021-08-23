#include "Traceback.hpp"

/*
TracebackResult *traceback(ScoringMatrix matrix, const char *a, const char *b) {
    char* a = calloc(sizeof(char), strlen(b));
}*/


std::vector<TracebackResult> traceback(ScoringMatrix matrix, const char *a, const char *b) {
    std::vector<std::pair<int,int>> mids; // store indices of maximum elementh
    mids.push_back(std::make_pair(0,0)); // init
    // find maximum value in scoring matrix
    for(int i=0;i<matrix.height;++i) {
        for(int j=0;j<matrix.width;++j) {
            if(matrix.scoring[(matrix.width * i + j)] > matrix.scoring[(matrix.width * mids[0].first + mids[0].second)]) {
                mids.clear();
                mids.push_back(std::make_pair(i,j));
            } else {
                if(matrix.scoring[(matrix.width * i + j)] == matrix.scoring[(matrix.width * mids[0].first + mids[0].second)]) {
                    mids.push_back(std::make_pair(i,j));
                }
            }
        }
    }

    std::vector<TracebackResult> res;
    for(auto &max : mids) {
        TracebackResult trc{"","",0,0,0,0.0,0.0};
        
        int x = max.first; // rows
        int y = max.second; // cols 
        int counter = 0;
        int matches = 0;

        std::string aAlign;
        std::string bAlign;

        while(matrix.scoring[matrix.width * x + y] > 0) {
            std::bitset<3> direction = matrix.traceback[matrix.width * x + y];
            if(direction.count() > 0) {
                if(direction[0] == 1) {
                    aAlign += a[x-1];
                    bAlign += b[y-1];
                    x--;
                    y--;
                    if(x > 0 && y > 0) {
                        if(matrix.scoring[matrix.width * x + y]-1 == matrix.scoring[matrix.width * (x-1) + (y-1)]) {
                            matches++;
                        }
                    }
                } else {
                    if(direction[1] == 1) {
                        aAlign += '-';
                        bAlign += b[y-1];
                        y--;
                    } else {
                        if(direction[2] == 1) {
                            aAlign += a[x-1];
                            bAlign += '-';
                            x--;
                        }
                    }
                }
            } else { // end reached
                break;
            }
            counter++;
        }

        std::reverse(aAlign.begin(),aAlign.end());
        std::reverse(bAlign.begin(),bAlign.end());
       
        // write back to object
        trc.a = aAlign;
        trc.b = bAlign;
        trc.length = aAlign.size();
        trc.matches = matches;
        trc.score = matrix.scoring[matrix.width * mids[0].first + mids[0].second];
        trc.cmpl = (double)trc.matches / trc.length;
        trc.ratio = (double)trc.matches / (strlen(a) + strlen(b));
        res.push_back(trc);
    }
    return res;
}




//
void free_traceback(TracebackResult trc) {
    
}



