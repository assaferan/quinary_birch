#ifndef __BIRCH_UTIL_H_
#define __BIRCH_UTIL_H_
// Birch-specific namespace for helper functions

#include "birch.h"

namespace birch_util
{
  template<typename From, typename To>
  To convert_Integer(const From& x);

  int popcnt(Z64 x);
}

#endif // __BIRCH_UTIL_H_
