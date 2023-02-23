#ifndef HELPER_HPP
#define HELPER_HPP

#include <boost/filesystem.hpp>
#include <iostream>
#include <iomanip>
#include <chrono>

namespace fs = boost::filesystem;

// filesystem manipulation
namespace helper {
    void createFolder(fs::path path);
    void deleteFolder(fs::path path); // deletes folder if its exists
    
    void splitFileN(fs::path source, fs::path targetPath, int sep);
    int lineCount(fs::path file, std::string command); // count lines of file
    std::vector<fs::path> listFiles(fs::path path);
    void concatFiles(fs::path target, std::vector<fs::path> files);

    std::string getTime(); // reports the current time


}

#endif // DATA_HPP
