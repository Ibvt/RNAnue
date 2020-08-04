#ifndef ClOSING_HPP
#define CLOSING_HPP

#include <vector>

class Closing {
    private:
        std::vector<std::string> quotes;
    public:
        Closing();
        std::vector<std::string> retrieveQuotes();
        void show(std::ostream& _str);
};

#endif // CLOSING_HPP
