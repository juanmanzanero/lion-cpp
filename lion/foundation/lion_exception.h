#ifndef LION_EXCEPTION_H
#define LION_EXCEPTION_H
#include <exception>
#include <string>
#include <iostream>

class lion_exception : public std::runtime_error
{
 public:
    using base_type = std::runtime_error;
    using base_type::base_type;
};

#endif
