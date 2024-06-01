#include "Stats.hpp"

Stats::Stats() : statsCounter(0), stats({}) {}

void Stats::setReadsCount(std::string condition, int increment) {
    {
        std::lock_guard<std::mutex> lock(statsMutex);
        stats[condition].readsCount += increment;
    }
}
void Stats::setAlignedCount(std::string condition, int increment) {
    {
        std::lock_guard<std::mutex> lock(statsMutex);
        stats[condition].alignedCount += increment;
    }
}
void Stats::setSplitsCount(std::string condition, int increment) {
    {
        std::lock_guard<std::mutex> lock(statsMutex);
        stats[condition].splitsCount += increment;
    }
}
void Stats::setMultSplitsCount(std::string condition, int increment) {
    {
        std::lock_guard<std::mutex> lock(statsMutex);
        stats[condition].multSplitsCount += increment;
    }
}
void Stats::setNSurvivedCount(std::string condition, int nsurvivedcount) {
    {
        std::lock_guard<std::mutex> lock(statsMutex);
        stats[condition].nSurvivedCount = nsurvivedcount;
    }
}