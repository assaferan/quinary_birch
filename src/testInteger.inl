#include <cassert>

template<typename R>
testInteger<R>::testInteger()
{
  testGcd();
  testXGcd();
}

template<typename R>
void testInteger<R>::testXGcd(const R & a, const R & b, const R & d,
			       const R & s, const R & t)
{
  Integer<R> a_int = a;
  Integer<R> b_int = b;

  typename EuclideanDomain< Integer<R> >::XGcdRes d_s_t = a_int.xgcd(b_int);
  
  Integer<R> d_int = std::get<0>(d_s_t);
  Integer<R> s_int = std::get<1>(d_s_t);
  Integer<R> t_int = std::get<2>(d_s_t);

  assert((d_int == d) && (s_int == s) && (t_int == t));

  return;
}

template<typename R>
void testInteger<R>::testGcd(const R & a, const R & b, const R & d)
{
  Integer<R> a_int = a;
  Integer<R> b_int = b;
  
  Integer<R> d_int = a_int.gcd(b_int);

  assert(d_int == d);

  return;
}

template<typename R>
void testInteger<R>::testXGcd()
{
  // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
  testXGcd(240, 46, 2, -9, 47);
  // https://en.wikipedia.org/wiki/Euclidean_algorithm
  testXGcd(252, 105, 21, -2, 5);
  testXGcd(1071, 462, 21, -3, 7);

  // test negative numbers
  testXGcd(-240, 46, 2, 9, 47);
  testXGcd(-240, -46, 2, 9, -47);
  
  return;
}

template<typename R>
void testInteger<R>::testGcd()
{
  // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
  testGcd(240, 46, 2);
  // https://en.wikipedia.org/wiki/Euclidean_algorithm
  testGcd(252, 105, 21);
  testGcd(1071, 462, 21);

  // test negative numbers
  testGcd(-240, 46, 2);
  testGcd(-240, -46, 2);

  return;
}
