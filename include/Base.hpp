#ifndef BASE_HPP
#define BASE_HPP

#include "Data.hpp"
#include "SeqRickshaw.hpp"

// namespace

class Base {
    private: 
        po::variables_map params;
        Data data;

    public:
        Base(po::variables_map params, std::string _subcall);

        template <typename Callable>
        void testBind(Callable f);
};

#endif // BASE_HPP
