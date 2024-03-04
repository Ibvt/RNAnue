//
// Created by Richard Albin Schaefer on 1/22/24.
//

#include <iostream>
#include <cstdlib>

#include "Closing.hpp"

Closing::Closing() : quotes(retrieveQuotes()) {}
Closing::~Closing() {}

std::vector<std::string> Closing::retrieveQuotes() {
    std::vector<std::string> quotes;
    quotes.insert(quotes.end(),("\"Hoppala\" - Cristobal"));
    quotes.insert(quotes.end(),("\"Neeeeinnnnnn\" - Biggi"));
    quotes.insert(quotes.end(),("\"Das ist nicht das Thema\" - Adruino "));
    quotes.insert(quotes.end(),("\"Haha lustig, genauso wie deine Kinder im Puff sehen\" - Biggi"));
    quotes.insert(quotes.end(),("\"Moment mal...\" - Randy"));
    quotes.insert(quotes.end(),("\"Wie du gehst schon?\" - Cristobal"));
    quotes.insert(quotes.end(),("\"Ist OK, OK?\" - Mr T."));
    quotes.insert(quotes.end(),("\"Number crunching, eye candy!\" - Adruino"));
    quotes.insert(quotes.end(),("\"Ich schiess mir...!\" - Anyone"));
    quotes.insert(quotes.end(),("\"Ich versteh das einfach nicht...!\" - Rich"));
    quotes.insert(quotes.end(),("\"Merci\" - Cristobal"));
    quotes.insert(quotes.end(),("\"Salut\" - Cristobal"));
    quotes.insert(quotes.end(),("\"Alla Jut\" - Cristobal"));
    quotes.insert(quotes.end(),("\"Make Baden great again\" - Cristobal"));
    quotes.insert(quotes.end(),("\"Mhm mhm mhm mhm\" - Cristobal"));
    quotes.insert(quotes.end(),("\"Ahhhhhhhh\" - Cristobal"));
    quotes.insert(quotes.end(),("\"Da isser ja (endlich)\" - Cristobal"));
    quotes.insert(quotes.end(),("\"Das ist kein Kaesekuchen, der ist kaesekuchenartig!\" - Cristobal"));
    quotes.insert(quotes.end(),("\"Come on!\" - Cristobal"));
    quotes.insert(quotes.end(),("\"Echt?!\" - Biggi"));
    quotes.insert(quotes.end(),("\"quasi gesehen\" - Biggi"));
    quotes.insert(quotes.end(),("\"Du Nazi\" - Biggi"));
    quotes.insert(quotes.end(),("\"Was stimmt mit dir nicht?\" - Biggi"));
    quotes.insert(quotes.end(),("\"Junger Vatter\" - Biggi"));
    quotes.insert(quotes.end(),("\"Ja voll!\" - Biggi"));
    quotes.insert(quotes.end(),("\"So dumm!\" - Biggi"));
    quotes.insert(quotes.end(),("\"Ja gut aehh...\" - Randy"));
    quotes.insert(quotes.end(),("\"In der frueh\" - Randy"));
    quotes.insert(quotes.end(),("\"*ruelps* ooopps...\" - Randy"));
    quotes.insert(quotes.end(),("\"oh leck ey\" - Randy"));
    quotes.insert(quotes.end(),("\"Ich mach das in Python\" - Randy"));
    quotes.insert(quotes.end(),("\"Flachzange\" - Randy"));
    quotes.insert(quotes.end(),("\"Das musst du wissen\" - Randy"));
    quotes.insert(quotes.end(),("\"Massephase\" - Randy"));
    quotes.insert(quotes.end(),("\"Spinnst du!\" - Randy"));
    quotes.insert(quotes.end(),("\"Hechtschenkel..?\" - Randy"));

    return quotes;
}

void Closing::printQuote(std::ostream& _str) {
    time_t t; //

    // intialize random number generator
    srand((unsigned) time(&t));

    // create random number and access quotes array
    int numberOfQuotes = this->quotes.size();
    int randomNumber = rand() % numberOfQuotes;
    // print random quote
    _str << this->quotes[randomNumber] << std::endl;
}






