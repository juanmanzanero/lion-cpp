#ifndef __TYPES_H__
#define __TYPES_H__

#include <cstddef>

using scalar = double;

#ifndef SKIP_TIMESERIES
using timeseries = double;  // For future integration with autodiff?
#endif 

using size_t = std::size_t;

#endif
