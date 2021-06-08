#include <cassert>

template<typename R>
testInteger<R>::testInteger()
{
  assert(testXGcd());
}

template<typename R>
bool testInteger<R>::testXGcd(const R & a, const R & b, const R & d,
			       const R & s, const R & t)
{
  Integer<R> a_int = a;
  Integer<R> b_int = b;

  typename EuclideanDomain< Integer<R> >::XGcdRes d_s_t = a_int.xgcd(b_int);
  
  Integer<R> d_int = std::get<0>(d_s_t);
  Integer<R> s_int = std::get<1>(d_s_t);
  Integer<R> t_int = std::get<2>(d_s_t);

  return ((d_int == d) && (s_int == s) && (t_int == t));
}

template<typename R>
bool testInteger<R>::testGcd(const R & a, const R & b, const R & d)
{
  Integer<R> a_int = a;
  Integer<R> b_int = b;
  
  Integer<R> d_int = a.gcd(b);

  return (d_int == d);
}

template<typename R>
bool testInteger<R>::testXGcd()
{
  // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
  return testXGcd(240, 46, 2, -9, -47);
  // https://en.wikipedia.org/wiki/Euclidean_algorithm
  return testXGcd(252, 105, 21, -2, 5);
  return testXGcd(1071, 462, 21, -3, 7);
}

template<typename R>
bool testInteger<R>::testGcd()
{
  // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
  return testGcd(240, 46, 2);
  // https://en.wikipedia.org/wiki/Euclidean_algorithm
  return testGcd(252, 105, 21);
  return testGcd(1071, 462, 21);
}
