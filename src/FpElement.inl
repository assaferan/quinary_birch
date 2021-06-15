template<typename R, typename S>
inline FpElement<R,S> FpElement<R,S>::operator-(void) const {
  assert(_GF != NULL);
  return FpElement(_GF, _GF->neg(_val));
}

template<typename R, typename S>
inline FpElement<R,S> FpElement<R,S>::operator+(const FpElement<R,S> &other) const
{
  assert(_GF != NULL);
  return FpElement(_GF, _GF->add(this->_val, other._val));
}

template<typename R, typename S>
inline FpElement<R,S> FpElement<R,S>::operator-(const FpElement<R,S> &other) const
{
  assert(_GF != NULL);
  return FpElement(_GF, _GF->sub(this->_val, other._val));
}

template<typename R, typename S>
inline FpElement<R,S> FpElement<R,S>::operator*(const FpElement<R,S> &other) const
{
  assert(_GF != NULL);
  return FpElement(_GF, _GF->mul(this->_val, other._val));
}

template<typename R, typename S>
inline FpElement<R,S> FpElement<R,S>::operator/(const FpElement<R,S> &other) const
{
  assert((_GF != NULL) && (other != 0));
  R elt = _GF->mod(other._val).lift();
  return FpElement(_GF, _GF->mul(this->_val, _GF->inverse(elt)));
}

template<typename R, typename S>
inline FpElement<R,S> FpElement<R,S>::inverse(void) const
{
  assert((_GF != NULL) && ((*this) != 0));
  R elt = _GF->mod(this->_val).lift();
  return FpElement(_GF, _GF->inverse(elt));
}

template<typename R, typename S>
inline FpElement<R,S> & FpElement<R,S>::operator+=(const FpElement<R,S> &other)
{
  assert(_GF != NULL);
  _val = _GF->add(this->_val, other._val); return (*this);
}

template<typename R, typename S>
inline FpElement<R,S> & FpElement<R,S>::operator-=(const FpElement<R,S> &other)
{
  assert(_GF != NULL);
  _val = _GF->sub(this->_val, other._val); return (*this);
}

template<typename R, typename S>
inline FpElement<R,S> & FpElement<R,S>::operator*=(const FpElement<R,S> &other)
{
  assert(_GF != NULL);
  _val = _GF->mul(this->_val, other._val); return (*this);
}

template<typename R, typename S>
inline FpElement<R,S> & FpElement<R,S>::operator/=(const FpElement<R,S> &other)
{
  assert((_GF != NULL) && (other != 0));
  R elt = _GF->mod(other._val).lift();
  _val = _GF->mul(this->_val, _GF->inverse(elt)); return (*this);
}

template<typename R, typename S>
inline int FpElement<R,S>::legendre(void) const
{
  return this->_GF->legendre(this->_val);
}

template<typename R, typename S>
inline FpElement<R,S> FpElement<R,S>::sqrt(void) const
{
  FpElement<R,S> a = *this;
  
  if (a.isOne()) return a;
  if (a.isZero()) return a;
  if (this->legendre() != 1) return this->_GF->zero();

  R p = this->_GF->prime();
  R q = p-1;
  R s = 0;
  while(q % 2 == 0) { q >>= 1; s++; }

  if (s == 1) return a^((p+1)/4);

  FpElement<R,S> z = this->_GF->random();
  while (z == 0 || z.legendre() == 1) {
    z = this->_GF->random();
  }

  int m = s;
  FpElement<R,S> c = z^q;
  FpElement<R,S> r = a^((q+1)/2);
  FpElement<R,S> t = a^q;

  while (1)
    {
      if (t.isOne()) return r;
      int i = 0;
      FpElement<R,S> t1 = t;
      while (!t1.isOne())
      {
	  t1 *= t1;
	  i++;
      }

      int e = 1;
      for (int j=0; j<m-i-1; j++) e <<= 1;
      FpElement<R,S> b = c^e;
      r *= b;
      c = b*b;
      t *= c;
      m = i;
    }

  return this->_GF->zero();
}

// assignment and conversion
template<typename R, typename S>
inline FpElement<R,S> & FpElement<R,S>::operator=(const FpElement<R,S> &other) 
{
  if (this != (&other)) {
    this->_GF = other._GF;
    this->_val = other._val;
  }
  return (*this);
}

//boolean
template<typename R, typename S>
inline bool FpElement<R,S>::operator==(const FpElement<R,S> &other) const {
  assert( (_GF != 0) && (other._GF != 0) ); 
  if (this->_GF->prime() != other._GF->prime()) return false;
  return ((this->_val - other._val) % (this->_GF->prime()) == 0);
}

template<typename R, typename S>
inline bool FpElement<R,S>::operator==(const R &other) const {
  assert(_GF != 0);
  return ((this->_val - other) % (this->_GF->prime()) == 0);
}

template<typename R, typename S>
inline bool FpElement<R,S>::operator!=(const R &other) const {
  assert(_GF != 0);
  return ((this->_val - other) % (this->_GF->prime()) != 0);
}

template<typename R, typename S>
inline bool FpElement<R,S>::isSquare(void) const {
  assert(_GF != 0);
  return ((this->_GF->legendre(this->_val)) >= 0);
}
