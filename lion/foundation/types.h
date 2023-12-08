#ifndef __TYPES_H__
#define __TYPES_H__

#include <cstddef>
#include <algorithm>

using scalar = double;

using std::size_t;
using std::max;
using std::min;


#ifdef _MSC_VER
#define windows_local_constexpr_auto static constexpr auto
#else
#define windows_local_constexpr_auto constexpr auto
#endif

#endif
