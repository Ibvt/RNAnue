#include "Base.hpp"

Base::Base(po::variables_map _params, std::string _subcall):
    params(_params),
    data(_params) {
        std::cout << "create Base" << std::endl;

        if(_params["subcall"].as<std::string>() == "preproc") {
            std::cout << "preproc called within Base" << std::endl;

            data.preproc();

            // create preprocessing instance (e.g., SeqRickshaw)
//            SeqRickshaw rsh;
 //           data.bla(std::bind(&SeqRickshaw::start, rsh));







//            std::cout << "call retrieveGroupsPath again" << std::endl;
 //           data.retrieveGroupsPath(fs::path(_params["ctrls"].as<std::string>()), fs::path(_params["trtms"].as<std::string>()));

//            data.callInAndOut(std::bind(&SeqRickshaw::start, rsh));


//            testBind(std::bind(&SeqRickshaw::start, rsh));
            // 
            //auto bla = std::bind(&SeqRickshaw::start, rsh, std::placeholders::_1);
            //pt::ptree something;
            //

 //           data.callInAndOut();

//            testBind(std::bind(&SeqRickshaw::start, rsh, std::placeholders::_1));
        }
    }

//
template <typename Callable>
void Base::testBind(Callable f) {
    std::cout << "asda" << std::endl;

    pt::ptree bla;

    f();
}

