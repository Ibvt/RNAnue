#ifndef RNANUE_DETECT_HPP
#define RNANUE_DETECT_HPP

// Standard
#include <iostream>
#include <bitset>
#include <stdlib.h>
#include <thread>
#include <mutex>
#include <future>
#include <vector>

// Boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

// HTSlib
#include <htslib/sam.h>
#include <htslib/hts.h>

// Vienna
//#include <ViennaRNA/fold.h>

// SeqAn3
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>

// Class
#include "IBPTree.hpp"
#include "Stats.hpp"

using types = seqan3::type_list<
        std::string,
        seqan3::sam_flag,
        std::optional<int32_t>,
        std::optional<int32_t>,
        std::optional<uint8_t>,
        dtp::CigarVector,
        dtp::DNASpan,
        seqan3::sam_tag_dictionary>;

using fields = seqan3::fields<
        seqan3::field::id,
        seqan3::field::flag,
        seqan3::field::ref_id,
        seqan3::field::ref_offset,
        seqan3::field::mapq,
        seqan3::field::cigar,
        seqan3::field::seq,
        seqan3::field::tags>;

using types2 = seqan3::type_list<
        std::vector<seqan3::dna5>,
        std::string,
        std::vector<seqan3::cigar>>;
using fields2 = seqan3::fields<
        seqan3::field::seq,
        seqan3::field::id,
        seqan3::field::cigar>;
using sam_record_type = seqan3::sam_record<types2, fields2>;




// define tags
using seqan3::operator""_tag;
template <> struct seqan3::sam_tag_type<"XH"_tag> { using type = int32_t; }; // number of splits (SAM record)
template <> struct seqan3::sam_tag_type<"XX"_tag> { using type = int32_t; }; // begin of split
template <> struct seqan3::sam_tag_type<"XY"_tag> { using type = int32_t; }; // end of split
template <> struct seqan3::sam_tag_type<"XJ"_tag> { using type = int32_t; }; // number of splits (whole read)
template <> struct seqan3::sam_tag_type<"XN"_tag> { using type = int32_t; };

// cigar



using SAMrecord = seqan3::sam_record<types, fields>;
using namespace seqan3::literals;
using seqan3::get;

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

class ThreadPool {
public:
    ThreadPool(size_t num_threads);
    ~ThreadPool();

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) -> std::future<std::invoke_result_t<F, Args...>>;

private:
    // Worker threads
    std::vector<std::thread> workers;
    // Task queue
    std::queue<std::function<void()>> tasks;

    // Synchronization
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};


class SplitReadCalling {
    public:
        SplitReadCalling(po::variables_map params);
        ~SplitReadCalling();

        void start(pt::ptree sample, pt::ptree condition);
        void sort(const std::string& inputFile, const std::string& outputFile);
        void iterate(std::string& matched, std::string& single, std::string& splits, std::string& multsplits);
        template <typename T>
        void process(std::vector<T>& readrecords, auto& singleOut, auto& splitsOut, auto& multsplitsOut,
                     auto& singleOutMutex, auto& splitsOutMutex, auto& multsplitsOutMutex);
        void decide(std::map<int, std::vector<SAMrecord>>& putative, auto& splitsOut, auto& multsplitsOut,
                    auto& splitsOutMutex, auto& multsplitsOutMutex);
        bool filter(auto& sequence, uint32_t cigarmatch, std::string chrom, dtp::Interval interval);
        void storeSegments(auto& splitrecord, std::optional<int32_t>& refOffset, dtp::CigarVector& cigar,
                           dtp::DNAVector& seq, seqan3::sam_tag_dictionary &tags, std::vector<SAMrecord>& curated);

    private:
        po::variables_map params;
        IBPTree features;
        //Stats stats;
        std::shared_ptr<Stats> stats;
        std::string condition; // stores the current condition
        std::deque<std::string> refIds; // stores the reference ids
};

#endif //RNANUE_DETECT_HPP
