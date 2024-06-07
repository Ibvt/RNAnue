#include "Utility.hpp"

// create folder
void helper::createDir(fs::path path, std::ostream& out) {
    if(fs::exists(path)) {
        out << getTime() << "The directory " << path << " already exists\n";
    } else {
        fs::create_directory(path);
        out << getTime() << "The directory " << path << " was created\n";
    }
}

void helper::createOutDir(fs::path path, std::ostream& out) {
    std::cout << getTime() << "Create directory to store the results (specified via --outdir)\n";
    createDir(path, out);
}

void helper::createTmpDir(fs::path path) {
    std::cout << helper::getTime() << "Create temporary directory " << path << "\n";
    deleteDir(path);
    fs::create_directory(path);
}

// lists the files in a directory
dtp::PathVector helper::listDirFiles(fs::path& path) {
    dtp::PathVector content;
    for(const auto& entry : fs::directory_iterator(path)) {

        if(entry.path().extension() != ".DS_Store") { // ignore tmp files in macos
            content.push_back(entry.path());
        }
    }

    // Sort paths using a lambda for case-insensitive comparison
    std::sort(content.begin(), content.end(), [](const fs::path& a, const fs::path& b) {
        auto cmp = [](unsigned char c1, unsigned char c2) {
            return std::tolower(c1) < std::tolower(c2);
        };
        return std::lexicographical_compare(a.string().begin(), a.string().end(),
                                            b.string().begin(), b.string().end(), cmp);
    });

    return content;
}

dtp::PathVector helper::filterDirFiles(dtp::PathVector& pathvec, std::string subcall) {
    dtp::PathVector filtered;
    for(auto& file : pathvec) {
        // in detect mode, only bam files are considered
        if(subcall == "detect") {
            if(file.extension() == ".bed" || file.extension() == ".txt") {
                continue;
            }
        }
        filtered.push_back(file);
    }
    return filtered;
}

bool helper::caseInsensitivePathCompare(const fs::path& a, const fs::path& b) {
    return a.string() < b.string();
}

void helper::renameFiles(fs::path dir, std::string extension) {
    dtp::PathVector files = helper::listDirFiles(dir);
    for(unsigned i=0;i<files.size();++i) {
        fs::path newPath = dir / fs::path(std::to_string(i) + extension);
        fs::rename(files[i], newPath);
        files[i] = newPath;
    }
}

std::string helper::removeNonPrintable(const std::string str) {
    std::string newStr;
    std::copy_if(str.begin(), str.end(), std::back_inserter(newStr), [](unsigned char c) {
        return std::isprint(c);
    });
    return newStr;
}


dtp::PathVector helper::genOutPath(dtp::PathVector inPathVec, std::string dir) {
    dtp::PathVector outPathVec;
    fs::path outPath = inPathVec[0].parent_path().parent_path() / fs::path(dir);

    for(auto& file : inPathVec) {
        outPathVec.push_back(outPath / fs::path(file.filename()));
    }
    return outPathVec;
}

// switch path (up to level condition from input to Output
fs::path helper::replacePath(fs::path _replacement, fs::path _original ) {
    fs::path nPath = _replacement / _original.filename();
    return nPath;
}


// delete folder
void helper::deleteDir(fs::path path) {
    if (fs::exists(path)) {
        fs::remove_all(path);
    }
}

// adds suffix to filename
std::string helper::addSuffix(std::string filename, std::string suffix, std::vector<std::string> keywords) {
    int keyPos, tmpKeyPos = -1; // buffer the positions of the keywords
    int dotPos = filename.find("."); // determine position of the dot
    keyPos = dotPos;
    if(!keywords.empty()) {
        for(unsigned i=0;i<keywords.size();++i) {
            tmpKeyPos = filename.find(keywords[i]);
            if(tmpKeyPos != -1) { // key could be found
                keyPos = tmpKeyPos;
            }
        }
    }
    std::string newFile = filename.substr(0,keyPos);
    newFile += suffix;
    newFile += filename.substr(dotPos,filename.size());

    return newFile;
}



void helper::splitFileN(fs::path source, fs::path targetPath, int sep) {
    fs::path targetStem = targetPath / "F";
    std::string call = "awk 'NR%" + std::to_string(sep) + "==1{x=\"" + targetStem.string();
    call += "\"++i\"" + source.extension().string() +"\";}{print>x}' " + source.string();

    system(call.c_str());
}

/*
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
}*/

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

bool helper::withinRange(int a, int b, int range) {
    if(std::abs(a-b) <= range || std::abs(b-a) <= range) {
        return true;
    } else {
        return false;
    }
}

std::string helper::getTime() {
    // Get the current time and date
    auto now = std::chrono::system_clock::now();
    std::time_t current_time = std::chrono::system_clock::to_time_t(now);

    // Convert the time and date to a string in the format: YYYY-MM-DD
    std::stringstream time_stream;
    time_stream << std::put_time(std::localtime(&current_time), "[%Y-%m-%d %H:%M:%S]") << " ";
    std::string time_string = time_stream.str();
    return time_string;
}

// ############ SEQIO ############
bool seqIO::filterReads(auto& qual, int quality, int minlen) {
    // determine the length
    size_t len = std::ranges::size(qual);
    if(calcAvgPhred(qual) >= quality && len >= minlen) {
        return true;
    } else {
        return false;
    }
}


dtp::DNAVector seqIO::spanToVector(dtp::DNASpan span) {
    dtp::DNAVector vec = std::vector<seqan3::dna5>{}; // empty vector
    vec.insert(vec.end(), span.begin(), span.end());
    return vec;
}

dtp::DNAVector seqIO::determineAlphabet(dtp::DNAVector seq) {
    std::vector<seqan3::dna5> alphabet{};
    // find the alphabet
    for (unsigned i = 0; i < seq.size(); ++i) {
        if (std::find(alphabet.begin(), alphabet.end(), seq[i]) == alphabet.end()) {
            alphabet.push_back(seq[i]);
        }
    }
    return alphabet;
}

void seqIO::printStates(dtp::DNAVector seq, std::string state, dtp::Left left, dtp::Right right, int readPos) {
    std::cout << "######################################\n";
    std::cout << "State: " << state << " | ";
    std::cout << "Left: " << left << " | ";
    std::cout << "Right: " << right.first << " " << right.second << " | ";
    std::cout << "ReadPos: " << readPos << " | ";
    std::cout << "Sequence: ";
    for(auto& s : seq) {
        std::cout << seqan3::to_char(s);
    }
    std::cout << "\n";
}

void seqIO::printDNAVector(dtp::DNAVector seq, std::ostream& ofs) {
    for(auto& s : seq) {
        ofs << seqan3::to_char(s);
    }
    ofs << "\n";
}

void seqIO::printDNASpan(dtp::DNASpan span, std::ostream& ofs) {
    for(auto& s : span) {
        ofs << seqan3::to_char(s);
    }
    ofs << "\n";
}

void seqIO::printQualSpan(dtp::QualSpan span, std::ostream& ofs) {
    for(auto& s : span) {
        ofs << seqan3::to_char(s);
    }
    ofs << "\n";
}

std::string seqIO::removeNonATGC(std::string seq) {
    std::string newSeq;
    for(auto& s : seq) {
        if(s == 'A' || s == 'T' || s == 'G' || s == 'C') {
            newSeq += s;
        }
    }
    return newSeq;
}


