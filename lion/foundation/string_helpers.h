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

inline std::string& strrep_in_place(std::string& str,
                                    const std::string& from,
                                    const std::string& to)
{
    //
    // Replaces all occurrences of "from" with "to"
    // in the input std::string "str", in place.
    // Returns a reference to the input string.
    //

    if (!from.empty()) {
        std::size_t pos{ 0u };
        while ((pos = str.find(from, pos)) != std::string::npos) {
            str.replace(pos, from.length(), to);
            pos += to.length(); // in case "to" contains "from", like replacing "x" with "yx"
        }
    }

    return str;
}

inline std::string strrep(std::string str,
                          const std::string& from,
                          const std::string& to)
{
    //
    // Replaces all occurrences of "from" with "to"
    // in the input std::string "str".
    //

    strrep_in_place(str, from, to);
    return str;
}

#endif