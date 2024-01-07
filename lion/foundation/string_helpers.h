#ifndef LION_FOUNDATION_STRING_HELPERS_H
#define LION_FOUNDATION_STRING_HELPERS_H
#pragma once


#include <cctype>
#include <string>
#include <sstream>
#include <limits>
#include <type_traits>


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


inline std::string& strleft_unquote_in_place(std::string &str)
{
    //
    // Removes the leftmost quotation marks (") from an
    // std::string, in-place. Returns a reference to the
    // input string.
    //

    str.erase(str.begin(), std::find_if(str.begin(), str.end(),
        [](unsigned char ch)
        { return ch != '"'; }));
    return str;
}

inline std::string& strright_unquote_in_place(std::string &str)
{
    //
    // Removes the rightmost quotation marks (") from an
    // std::string, in-place. Returns a reference to the
    // input string.
    //

    str.erase(std::find_if(str.rbegin(), str.rend(),
        [](unsigned char ch)
        { return ch != '"'; }).base(), str.end());
    return str;
}

inline std::string& strunquote_in_place(std::string &str)
{
    //
    // Removes the leftmost and rightmost quotation marks (")
    // from an std::string, in-place. Returns a reference to
    // the input string.
    //

    strleft_unquote_in_place(str);
    strright_unquote_in_place(str);
    return str;
}

inline std::string strunquote(std::string str)
{
    //
    // Returns the input std::string without its
    // leftmost and rightmost quotation marks.
    //

    strunquote_in_place(str);
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


template<typename T,
         typename std::enable_if_t<std::is_arithmetic_v<T> >* = nullptr>
inline std::string num2str(const T &num)
{
    //
    // Does the "number -> std::string" conversion,
    // maintaining the precision of the number's type.
    //

    if constexpr (std::is_integral_v<T>) {
        // integral types don't have decimal places
        return std::to_string(num);
    }
    else {
        // floating-point types go with their max. decimal precision
        static_assert(std::is_floating_point_v<T>);

        std::ostringstream ss;
        ss.precision(std::numeric_limits<T>::max_digits10);
        ss << num;
        return ss.str();
    }
}

template<typename It>
inline std::string range2str(It first, It last,
                             const std::string &separator = ", ")
{
    //
    // Does the "[first, last) -> std::string" conversion,
    // maintaining the precision of the range's value_type,
    // with an optional custom separator.
    //

    if (first == last) {
        return std::string{};
    }
    else {
        std::ostringstream ss;
        ss.precision(std::numeric_limits<typename std::iterator_traits<It>::value_type>::max_digits10);
        ss << *first++;
        while (first != last) {
            ss << separator << *first++;
        }
        return ss.str();
    }
}

template<typename Vec>
inline std::string vec2str(const Vec &vec,
                           const std::string &separator = ", ")
{
    //
    // Does the "vector -> std::string" conversion,
    // maintaining the precision of the vector's value_type,
    // with an optional custom separator. This vector
    // may be an std::vector, an std::array, etc.
    //

    return range2str(vec.cbegin(), vec.cend(), separator);
}

#endif