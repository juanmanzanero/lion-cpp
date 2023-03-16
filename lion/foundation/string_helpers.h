#ifndef LION_FOUNDATION_STRING_HELPERS_H
#define LION_FOUNDATION_STRING_HELPERS_H
#pragma once


#include <cctype>
#include <string>


//
// Defines utility fuctions to deal with std::strings.
//


inline std::string strtolower(std::string str)
{
    //
    // Converts the input string to lowercase.
    //

    std::transform(str.cbegin(), str.cend(),
        str.begin(), [](char c)
                     {
                         return static_cast<char>(std::tolower(
                             static_cast<unsigned char>(c)));
                     });
    return str;
}


inline std::string strtoupper(std::string str)
{
    //
    // Converts the input string to uppercase.
    //

    std::transform(str.cbegin(), str.cend(),
        str.begin(), [](char c)
                     {
                         return static_cast<char>(std::toupper(
                             static_cast<unsigned char>(c)));
                     });
    return str;
}


inline std::string& strleft_trim_in_place(std::string &str)
{
    //
    // Left trims an std::string, in-place. Returns
    // a reference to the input string.
    //

    str.erase(str.begin(), std::find_if(str.begin(), str.end(),
        [](unsigned char ch)
        { return !std::isspace(ch); }));
    return str;
}


inline std::string& strright_trim_in_place(std::string &str)
{
    //
    // Right-trims an std::string, in-place. Returns
    // a reference to the input string.
    //

    str.erase(std::find_if(str.rbegin(), str.rend(),
        [](unsigned char ch)
        { return !std::isspace(ch); }).base(), str.end());
    return str;
}


inline std::string& strtrim_in_place(std::string &str)
{
    //
    // Trims an std::string, in-place. Returns
    // a reference to the input string.
    //

    strleft_trim_in_place(str);
    strright_trim_in_place(str);
    return str;
}


inline std::string strtrim(std::string str)
{
    //
    // Returns the trimmed version of the
    // input std::string.
    //

    strtrim_in_place(str);
    return str;
}

#endif
