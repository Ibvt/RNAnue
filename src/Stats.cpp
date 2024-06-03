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

void Stats::writeStats(fs::path outdir) {
    fs::path statsFile = outdir / "stats.tsv";
    std::ofstream statsStream(statsFile);
    statsStream << "Condition\tReads\tAligned\tSplits\tMultSplits\n";
    for (const auto& [condition, stat] : stats) {
        statsStream << condition << "\t" << stat.readsCount << "\t";
        statsStream << stat.alignedCount << "\t" << stat.splitsCount << "\t";
        statsStream << stat.multSplitsCount << "\n";
    }
}