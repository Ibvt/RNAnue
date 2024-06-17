#ifndef RNANUE_STATS_HPP
#define RNANUE_STATS_HPP

// Standard
#include <map>
#include <string>
#include <vector>
#include <mutex>

// Boost
#include <boost/filesystem.hpp>

// Class
#include "DataTypes.hpp"

namespace fs = boost::filesystem;

class Stats {
    public:
        Stats();
        Stats(std::string statsFile);

        // Move constructor / no needed because its unique
        Stats(Stats&& other) noexcept : stats(other.stats) {}

        // Move assignment operator
        Stats& operator=(Stats&& other) noexcept {
            if (this != &other) {
                std::lock_guard<std::mutex> lock(statsMutex);
                stats = other.stats;
            }
            return *this;
        }

        // Delete copy constructor and copy assignment operator
        Stats(const Stats&) = delete;
        Stats& operator=(const Stats&) = delete;

        // getter & setter
        void setReadsCount(std::string condition, int repl, int increment);
        void setAlignedCount(std::string condition, int repl, int increment);
        void setSplitsCount(std::string condition, int repl, int increment);
        void setMultSplitsCount(std::string condition, int repl, int increment);
        void setInteractionsCount(std::string condition, int repl, int increment);

        void reserveStats(std::string condition, int repl); // creates new entry for replicate

        // write stats back to file
        void writeStats(fs::path outdir, std::string subcall);

    private:
        dtp::StatsMap stats;
        int statsCounter; // counter for each library
        std::mutex statsMutex;
};

#endif //RNANUE_STATS_HPP
