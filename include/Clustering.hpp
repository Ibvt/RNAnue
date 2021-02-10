// boost
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

namespace po = boost::program_options;
namespace pt = boost::property_tree;

class Clustering {
    private:
        po::variables_map params;

	public:
		Clustering(po::variables_map params);
		Clustering();


        void start(pt::ptree sample);
};

