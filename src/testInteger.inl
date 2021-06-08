#include <cassert>

template<typename R>
testInteger<R>::testInteger()
{
  std::cerr << "testing gcd..." << std::endl;
  std::cerr << "testGcd() == " << testGcd() << std::endl;
  assert(testGcd() == true);
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
  
  Integer<R> d_int = a_int.gcd(b_int);

  return (d_int == d);
}

template<typename R>
bool testInteger<R>::testXGcd()
{
  bool ret = true;
  // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
  ret = ret && testXGcd(240, 46, 2, -9, -47);
  // https://en.wikipedia.org/wiki/Euclidean_algorithm
  ret = ret && testXGcd(252, 105, 21, -2, 5);
  ret = ret && testXGcd(1071, 462, 21, -3, 7);

  // test negative numbers
  ret = ret && testXGcd(-240, 46, 2, 9, -47);
  ret = ret && testXGcd(-240, -46, 2, 9, 47);
  
  return ret;
}

template<typename R>
bool testInteger<R>::testGcd()
{
  bool ret = true;
  // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
  ret = ret && testGcd(240, 46, 2);
  // https://en.wikipedia.org/wiki/Euclidean_algorithm
  ret = ret && testGcd(252, 105, 21);
  ret = ret && testGcd(1071, 462, 21);

  // test negative numbers
  ret = ret && testGcd(-240, 46, 2);
  ret = ret && testGcd(-240, -46, 2);

  ret = ret && false;
  return ret;
}
