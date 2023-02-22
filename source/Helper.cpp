#include "Helper.hpp"

// create folders - and delete
void helper::createFolder(fs::path path) {
    if(!fs::exists(path)) {
        fs::create_directory(path);
    }
}

void helper::deleteFolder(fs::path path) {
    if(fs::exists(path)) {
        fs::remove_all(path);
    }
}

// 
void helper::splitFileN(fs::path source, fs::path targetPath, int sep) {
    fs::path targetStem = targetPath / "F";
    std::string call = "awk 'NR%" + std::to_string(sep) + "==1{x=\"" + targetStem.string();
    call += "\"++i\"" + source.extension().string() +"\";}{print>x}' " + source.string();

    system(call.c_str());
}


int helper::lineCount(fs::path file, std::string command) {
    std::string entries = "cat " + file.string() + " " + command + " | wc -l";
    const char* call = entries.c_str();

    std::array<char, 128> buffer;

    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(call, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return std::stoi(result);
}


std::vector<fs::path> helper::listFiles(fs::path path) {
    std::vector<fs::path> files; // 
    copy(fs::directory_iterator(path), fs::directory_iterator(), back_inserter(files));
    sort(files.begin(),files.end()); // sort the content 
    return files;
}
    
void helper::concatFiles(fs::path target, std::vector<fs::path> files) {
    std::string call = "cat";
    for(unsigned i=0;i<files.size();++i) {
        call += " " + files[i].string();
    }
    call += " > " + target.string();
    if(files.size() > 0) {
        system(call.c_str());
    }
}

//
std::string helper::getTime() {
    // Get the current time and date
    auto now = std::chrono::system_clock::now();
    std::time_t current_time = std::chrono::system_clock::to_time_t(now);

    // Convert the time and date to a string in the format: YYYY-MM-DD
    std::stringstream time_stream;
    time_stream << std::put_time(std::localtime(&current_time), "[%Y-%m-%d %H:%M:%S]");
    std::string time_string = time_stream.str();
    return time_string;
}




