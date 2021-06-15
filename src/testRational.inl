#include <random>

template<typename R>
testRational<R>::testRational()
{
  R num1 = std::rand();
  R denom1 = std::rand();
  R num2 = std::rand();
  R denom2 = std::rand();

  Z64 e = std::rand();
  
  _testConstructor(num1, denom1);
  _testConstructor(num2, denom2);
  _testAdd(num1, denom1, num2, denom2);
  _testSub(num1, denom1, num2, denom2);
  _testMul(num1, denom1, num2, denom2);
  _testDiv(num1, denom1, num2, denom2);
  _testPow(num1, denom1, e);
  _testInverse(num1, denom1);
  _testZero(denom1);
  _testOne(num1);
  _testGcd(num1, denom1, num2, denom2);
  _testEuclid(num1, denom1, num2, denom2);
  
}

template<typename R>
inline void testRational<R>::_testConstructor(const R & num, const R & denom)
{
#ifdef DEBUG
  Rational<R> f(num, denom);
#endif
  assert(denom*f.num().num() == num*f.denom().num());
  return;
}

template<typename R>
inline void testRational<R>::_testAdd(const R & num1, const R & denom1,
			       const R & num2, const R & denom2)
{
#ifdef DEBUG
  Rational<R> f1(num1, denom1);
  Rational<R> f2(num2, denom2);

  Rational<R> f = f1 + f2;
#endif
  assert(f1.denom()*f2.denom()*f.num() ==
	 f.denom()*f2.denom()*f1.num() + f.denom()*f1.denom()*f2.num());

  return;
}

template<typename R>
inline void testRational<R>::_testSub(const R & num1, const R & denom1,
			       const R & num2, const R & denom2)
{
#ifdef DEBUG
  Rational<R> f1(num1, denom1);
  Rational<R> f2(num2, denom2);

  Rational<R> f = f1 - f2;
#endif
  assert(f1.denom()*f2.denom()*f.num() ==
	 f.denom()*f2.denom()*f1.num() - f.denom()*f1.denom()*f2.num());

  return;
}

template<typename R>
inline void testRational<R>::_testMul(const R & num1, const R & denom1,
			       const R & num2, const R & denom2)
{
#ifdef DEBUG
  Rational<R> f1(num1, denom1);
  Rational<R> f2(num2, denom2);

  Rational<R> f = f1 * f2;
#endif
  assert(f1.denom()*f2.denom()*f.num() ==  f.denom()*f1.num()*f2.num());

  return;
}

template<typename R>
inline void testRational<R>::_testDiv(const R & num1, const R & denom1,
			       const R & num2, const R & denom2)
{
#ifdef DEBUG
  Rational<R> f1(num1, denom1);
  Rational<R> f2(num2, denom2);

  Rational<R> f = f1 / f2;
#endif
  assert(f1.denom()*f2.num()*f.num() ==  f.denom()*f1.num()*f2.denom());

  return;
}

template<typename R>  
inline void testRational<R>::_testPow(const R & num1, const R & denom1, const Z64 & e)
{
#ifdef DEBUG
  Rational<R> f(num1, denom1);
  Rational<R> f_e = f^e;
#endif
  assert(f_e.num()*(f.denom()^e) == f_e.denom()*(f.num()^e));

  return;
}

template<typename R>  
inline void testRational<R>::_testInverse(const R & num, const R & denom)
{
#ifdef DEBUG
  Rational<R> f(num, denom);
  Rational<R> f_inv = f.inverse();
#endif
  assert((f*f_inv).isOne() && (f_inv*f).isOne());

  return;
}

template<typename R>
inline void testRational<R>::_testZero(const R & denom)
{
  Rational<R> f(0, denom);
  assert(f.isZero());
  f.makeZero();
  assert(f.isZero());

  return;
}

template<typename R>
inline void testRational<R>::_testOne(const R & val)
{
  Rational<R> f(val, val);
  assert(f.isOne());
  f.makeOne();
  assert(f.isOne());

  return;
}

template<typename R>
inline void testRational<R>::_testGcd(const R & num1, const R & denom1,
			       const R & num2, const R & denom2)
{
#ifdef DEBUG
  Rational<R> f1(num1, denom1);
  Rational<R> f2(num2, denom2);

  Rational<R> d = f1.gcd(f2);
#endif
  
  assert(d.isOne());

  return;
}

template<typename R>
inline void testRational<R>::_testEuclid(const R & num1, const R & denom1,
				  const R & num2, const R & denom2)
{
#ifdef DEBUG
  Rational<R> f1(num1, denom1);
  Rational<R> f2(num2, denom2);

  std::pair<Rational<R>, Rational<R> > q_r = f1.euclideanDivision(f2);

  assert(q_r.second.isZero());

  Rational<R> r = f1 % f2;

  assert(r.isZero());
#endif
  
  return;
}
  
