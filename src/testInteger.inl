#include <cassert>

template<typename R>
testInteger<R>::testInteger()
{
  _testGcd();
  _testXGcd();
}

template<typename R>
inline void testInteger<R>::_testXGcd(const R & a, const R & b, const R & d,
			       const R & s, const R & t)
{
#ifdef DEBUG
  Integer<R> a_int = a;
  Integer<R> b_int = b;

  typename EuclideanDomainElement< Integer<R>, IntegerRing<R> >::XGcdRes
    d_s_t = a_int.xgcd(b_int);
  
  Integer<R> d_int = std::get<0>(d_s_t);
  Integer<R> s_int = std::get<1>(d_s_t);
  Integer<R> t_int = std::get<2>(d_s_t);
#endif
  
  assert((d_int == d) && (s_int == s) && (t_int == t));

  return;
}

template<typename R>
inline void testInteger<R>::_testGcd(const R & a, const R & b, const R & d)
{
  
#ifdef DEBUG
  Integer<R> a_int = a;
  Integer<R> b_int = b;
  
  Integer<R> d_int = a_int.gcd(b_int);
#endif
  
  assert(d_int == d);

  return;
}

template<typename R>
inline void testInteger<R>::_testXGcd(void)
{
  // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
  _testXGcd(240, 46, 2, -9, 47);
  // https://en.wikipedia.org/wiki/Euclidean_algorithm
  _testXGcd(252, 105, 21, -2, 5);
  _testXGcd(1071, 462, 21, -3, 7);

  // test negative numbers
  _testXGcd(-240, 46, 2, 9, 47);
  _testXGcd(-240, -46, 2, 9, -47);
  
  return;
}

template<typename R>
inline void testInteger<R>::_testGcd(void)
{
  // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
  _testGcd(240, 46, 2);
  // https://en.wikipedia.org/wiki/Euclidean_algorithm
  _testGcd(252, 105, 21);
  _testGcd(1071, 462, 21);

  // test negative numbers
  _testGcd(-240, 46, 2);
  _testGcd(-240, -46, 2);

  return;
}
