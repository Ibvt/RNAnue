#include "Stats.hpp"

Stats::Stats() : statsCounter(0), stats({}) {}
Stats::Stats(std::string statsFile) {
    this->statsCounter = 0;
    this->stats = {};

    // read stats from file
    std::ifstream statsStream(statsFile);
    std::string line;
    if(std::getline(statsStream, line)) {} // skip header
    while(std::getline(statsStream, line)) {
        std::string token;
        std::vector<std::string> tokens;
        std::istringstream ss(line);
        while(getline(ss, token, '\t')) { tokens.push_back(token); }

        std::string condition = tokens[0];
        std::string readsCount = tokens[2];
        std::string alignedCount = tokens[3];
        std::string splitsCount = tokens[4];
        std::string multSplitsCount = tokens[5];

        dtp::StatsFields entry(std::stoi(readsCount), std::stoi(alignedCount),
                               std::stoi(splitsCount), std::stoi(multSplitsCount));

        // check if condition exists
        if(stats.find(condition) == stats.end()) {
            std::vector<dtp::StatsFields> tmp;
            tmp.push_back(entry);
            stats.insert({condition, tmp});
        } else {
            stats[condition].push_back(entry);
        }
    }
}

void Stats::setReadsCount(std::string condition, int repl, int increment) {
    {
        std::lock_guard<std::mutex> lock(statsMutex);
        stats[condition][repl].readsCount += increment;
    }
}
void Stats::setAlignedCount(std::string condition, int repl, int increment) {
    {
        std::lock_guard<std::mutex> lock(statsMutex);
        stats[condition][repl].alignedCount += increment;
    }
}
void Stats::setSplitsCount(std::string condition, int repl, int increment) {
    {
        std::lock_guard<std::mutex> lock(statsMutex);
        stats[condition][repl].splitsCount += increment;
    }
}
void Stats::setMultSplitsCount(std::string condition, int repl, int increment) {
    {
        std::lock_guard<std::mutex> lock(statsMutex);
        stats[condition][repl].multSplitsCount += increment;
    }
}

void Stats::setInteractionsCount(std::string condition, int repl, int increment) {
    {
        std::lock_guard<std::mutex> lock(statsMutex);
        stats[condition][repl].interactionsCount += increment;
    }
}

void Stats::reserveStats(std::string condition, int repl) {
    // check if condition exists
    if(stats.find(condition) == stats.end()) {
        {
            std::lock_guard<std::mutex> lock(statsMutex);
            std::vector<dtp::StatsFields> tmp;
            tmp.push_back(dtp::StatsFields());
            stats.insert({condition, tmp});
        }
    } else {
        if(stats[condition].size() <= repl) {
            {
                std::lock_guard<std::mutex> lock(statsMutex);
                stats[condition].push_back(dtp::StatsFields());
            }
        }
    }
}

void Stats::writeStats(fs::path outdir, std::string subcall) {
    fs::path statsFile = outdir / "stats.txt";
    std::ofstream statsStream(statsFile);
    statsStream << "Condition\tReplicate\tReads\tAligned\tSplits\tMultSplits";
    if(subcall == "analysis") { statsStream << "\tInteractions"; }
    statsStream << "\n";
    for (const auto& cond: stats) {
        for(int i=0; i<cond.second.size(); ++i) {
            std::stringstream oss;
            oss << std::setw(3) << std::setfill('0') << i;
            statsStream << cond.first << "\t" << oss.str();
            statsStream << "\t" << cond.second[i].readsCount;
            statsStream << "\t" << cond.second[i].alignedCount;
            statsStream << "\t" << cond.second[i].splitsCount;
            statsStream << "\t" << cond.second[i].multSplitsCount;
            if(subcall == "analysis") {
                statsStream << "\t" << cond.second[i].interactionsCount;
            }
            statsStream << "\n";
        }
    }
    statsStream.close();
}