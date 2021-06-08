#ifndef __TEST_INTEGERS_H_
#define __TEST_INTEGERS_H_

#include "Integer.h"

template<typename R>
class testInteger
{
public:
  testInteger();
  
protected:

  static bool testXGcd();
  static bool testGcd();
  
  static bool testXGcd(const R & a, const R & b, const R & d,
		       const R & s, const R & t);
  static bool testGcd(const R & a, const R & b, const R & d);
  
};

#include "testInteger.inl"

#endif //  __TEST_INTEGERS_H_
