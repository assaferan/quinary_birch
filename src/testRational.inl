#include <random>

template<typename R>
testRational<R>::testRational()
{
  R num1 = std::rand();
  R denom1 = std::rand();
  R num2 = std::rand();
  R denom2 = std::rand();

  Z64 e = std::rand();
  
  testConstructor(num1, denom1);
  testConstructor(num2, denom2);
  testAdd(num1, denom1, num2, denom2);
  testSub(num1, denom1, num2, denom2);
  testMul(num1, denom1, num2, denom2);
  testDiv(num1, denom1, num2, denom2);
  testPow(num1, denom1, e);
  testInverse(num1, denom1);
  testZero(denom1);
  testOne(num1);
  testGcd(num1, denom1, num2, denom2);
  testEuclid(num1, denom1, num2, denom2);
  
}

template<typename R>
void testRational<R>::testConstructor(const R & num, const R & denom)
{
#ifdef DEBUG
  Rational<R> f(num, denom);
#endif
  assert(denom*f.num().num() == num*f.denom().num());
  return;
}

template<typename R>
void testRational<R>::testAdd(const R & num1, const R & denom1,
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
void testRational<R>::testSub(const R & num1, const R & denom1,
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
void testRational<R>::testMul(const R & num1, const R & denom1,
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
void testRational<R>::testDiv(const R & num1, const R & denom1,
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
void testRational<R>::testPow(const R & num1, const R & denom1, const Z64 & e)
{
#ifdef DEBUG
  Rational<R> f(num1, denom1);
  Rational<R> f_e = f^e;
#endif
  assert(f_e.num()*(f.denom()^e) == f_e.denom()*(f.num()^e));

  return;
}

template<typename R>  
void testRational<R>::testInverse(const R & num, const R & denom)
{
#ifdef DEBUG
  Rational<R> f(num, denom);
  Rational<R> f_inv = f.inverse();
#endif
  assert((f*f_inv).isOne() && (f_inv*f).isOne());

  return;
}

template<typename R>
void testRational<R>::testZero(const R & denom)
{
  Rational<R> f(0, denom);
  assert(f.isZero());
  f.makeZero();
  assert(f.isZero());

  return;
}

template<typename R>
void testRational<R>::testOne(const R & val)
{
  Rational<R> f(val, val);
  assert(f.isOne());
  f.makeOne();
  assert(f.isOne());

  return;
}

template<typename R>
void testRational<R>::testGcd(const R & num1, const R & denom1,
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
void testRational<R>::testEuclid(const R & num1, const R & denom1,
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
  
