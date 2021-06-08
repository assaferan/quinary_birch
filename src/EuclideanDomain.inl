

template<class Derived>
typename EuclideanDomain<Derived>::XGcdRes
EuclideanDomain<Derived>::xgcd(const Derived& b) const
{
  if (b.isZero()) {
    // We don't know whether there exists a default constructor,
    // but there is a copy constructor.
    Derived a(*static_cast<const Derived *>(this));
    Derived d = a;
    Derived s = d;
    Derived t = d;
    s.one();
    t.zero();
    return std::make_tuple<Derived>(std::move(d),s,t);
  }

  typename EuclideanDomain<Derived>::DivRes q_r = this->euclideanDivision(b);

  Derived q = q_r.first;
  Derived r = q_r.second;

  // !! TODO - eliminate recursion here 
  typename EuclideanDomain<Derived>::XGcdRes d_t_s = b.xgcd(r);
  Derived d = std::get<0>(d_t_s);
  Derived t = std::get<1>(d_t_s);
  Derived s = std::get<2>(d_t_s);
  t -= q * s;
  
  return std::make_tuple<Derived>(d,s,t);
}

// This could use xgcd, but this implementation is lightly quicker
template<class Derived>
Derived EuclideanDomain<Derived>::gcd(const Derived& b) const
{
  if (b.isZero())
    return (*this);

  typename EuclideanDomain<Derived>::DivRes q_r = this->euclideanDivision(b);

  Derived q = q_r.first;
  Derived r = q_r.second;
  
  // !! TODO - eliminate recursion here
  return b.gcd(r);
}
  
/**
 * Get the quotient of *this and b.
 *
 * @param b: the divisor
 * @return the quotient
 */
template<class Derived>
Derived EuclideanDomain<Derived>::quotient(const Derived& b) const
{
  typename EuclideanDomain<Derived>::DivRes q_r = this->euclideanDivision(b);
  return q_r.first;
}

/**
 * Get the remainder of *this and b.
 * @param b: the divisor
 * @return the remainder
 */
template<class Derived>
Derived EuclideanDomain<Derived>::remainder(const Derived& b) const
{
  typename EuclideanDomain<Derived>::DivRes q_r = this->euclideanDivision(b);
  return q_r.second;
}

/**
 * Exact division.
 *
 * @param d: the divisor.
 * @return the equotient.
 */
template<class Derived>
Derived EuclideanDomain<Derived>::operator/ (const Derived& d) const
{
  return this->quotient(d);
}

/**
 * Exact division assignment.
 *
 * @param d: the divisor.
 * @return a reference to this after assignment.
 */
template<class Derived>
Derived& EuclideanDomain<Derived>::operator/= (const Derived& d)
{
  return ((*this) = (*this) / d);
}

/**
 * Get the remainder of *this and b;
 * @param b: the divisor
 * @return the remainder
 */
template<class Derived>
Derived EuclideanDomain<Derived>::operator%(const Derived& b) const
{
  return this->remainder(b);
}

/**
 * Assign *this to be the remainder of *this and b.
 * @param b: the divisor
 * @return this after assignment.
 */
template<class Derived>
Derived& EuclideanDomain<Derived>::operator%=(const Derived& b)
{
  return ((*this) = (*this) % b);
}
