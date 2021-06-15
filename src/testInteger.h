#ifndef __TEST_INTEGERS_H_
#define __TEST_INTEGERS_H_

#include "Integer.h"

template<typename R>
class testInteger
{
public:
  testInteger();
  
protected:

  static void _testXGcd(void);
  static void _testGcd(void);
  
  static void _testXGcd(const R & a, const R & b, const R & d,
		       const R & s, const R & t);
  static void _testGcd(const R & a, const R & b, const R & d);
  
};

#include "testInteger.inl"

#endif //  __TEST_INTEGERS_H_
