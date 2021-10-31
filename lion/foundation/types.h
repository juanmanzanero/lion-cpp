#ifndef __TYPES_H__
#define __TYPES_H__

#include <cstddef>

using scalar = double;

#ifndef SKIP_TIMESERIES
using timeseries = double;  // For future integration with autodiff?
#endif 

using std::size_t;
using std::max;
using std::min;

#endif
