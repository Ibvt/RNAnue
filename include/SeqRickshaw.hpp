#ifndef SEQRICKSHAW_H
#define SEQRICKSHAW_H

#include <bitset>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <algorithm>
#include <fstream>

// boost
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>

// seqan
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/io/sequence_file/format_fastq.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp> 
#include <seqan3/io/exception.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/slice.hpp>  

#include <seqan3/core/debug_stream.hpp>

#include <seqan3/std/ranges>

#include <seqan3/range/views/all.hpp>
#include <seqan3/std/ranges> // std::ranges::copy

#include <range/v3/all.hpp>
#include <range/v3/view/transform.hpp>


namespace pt = boost::property_tree;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

typedef std::tuple<std::string,int,std::pair<int,int>,std::size_t> State;
typedef std::vector<State> States;

typedef std::map<std::pair<int,char>,std::tuple<int,int,int,int>> LookupTable;

typedef std::map<std::pair<std::string,std::string>, LookupTable> Adapters;

//typedef std::map<std::pair<int,string::string>
//

using seqan3::operator""_dna5;
using seqan3::operator""_dna4;

class SeqRickshaw {
    private:
        int modus; // defines the modus of the trimming procedure (

        std::string adpt5prime;
        std::string adpt3prime;
        
        Adapters adpt5Table;
        Adapters adpt3Table;

        int phred; // the minimum 
        int minlen;

        std::string readtype;

        /* the lookup table - smart transition table 
         * (state,c) -> shift, state, readPos, match
         * */
//        std::map<std::pair<int,char>,std::tuple<int,int,int,int>> lookup;
            
    public:
        SeqRickshaw(po::variables_map _params);
        SeqRickshaw();

        std::map<std::pair<std::string,std::string>,LookupTable> calcLookupTable(std::string _type, std::string _path);
        LookupTable calcShift(auto _records);

        std::vector<char> determineAlphabet(auto _sequence);

       
        std::size_t calcReadPos(auto& sequence, std::size_t& left, std::pair<std::size_t,std::size_t>& right);



        void smallestShift(std::string pattern, std::string suffix, int left);

        int transition(std::string pattern, std::string suffix, int readPos, std::size_t& left, std::pair<std::size_t,std::size_t>& right);

        void merging(auto fwd, auto rev);


        std::string longestCommonSubstr(std::string forward, std::string reverse);




        // helper 
        int addState(States &states, State state, States::size_type &size);
        int nextReadPos(std::string state, int currReadPos);
        // finds all occurrences of substring in string
        void findAllOcc(std::vector<std::size_t>& fnd, std::string str, std::string substr);
        void writeLookupTable(std::ofstream &os);

        void preprocPattern();
        
        
        std::size_t boyermoore(auto& read, LookupTable tab, int patlen);
        std::pair<std::size_t,std::size_t> trimming(auto& read);
        void start(pt::ptree sample);

        bool filtering(auto &rec);

};

#endif
