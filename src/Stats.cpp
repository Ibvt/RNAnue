#include "Stats.hpp"

Stats::Stats() : statsCounter(0), stats({}) {
    // initialize stats counter
    statsCounter = 0;

}
/*
void Stats::setReadsCount(int readscount) { stats["readscount"] = readscount; }
void Stats::setAlignedCount(int alignedcount) { stats["alignedcount"] = alignedcount; }
void Stats::setSplitsCount(int splitscount) { stats["splitscount"] = splitscount; }
void Stats::setMultSplitsCount(int multsplitscount) { stats["msplitscount"] = multsplitscount; }
void Stats::setNSurvivedCount(int nsurvivedcount) { stats["nsurvivedcount"] = nsurvivedcount; }
 */