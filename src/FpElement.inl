template<typename R, typename S>
FpElement<R,S> FpElement<R,S>::operator-() const {
  assert(_GF != NULL);
  return FpElement(_GF, _GF->neg(_val));
}

template<typename R, typename S>
FpElement<R,S> FpElement<R,S>::operator+(const FpElement<R, S> &other) const
{
  assert(_GF != NULL);
  return FpElement(_GF, _GF->add(this->_val, other._val));
}

template<typename R, typename S>
FpElement<R,S> FpElement<R,S>::operator-(const FpElement<R, S> &other) const
{
  assert(_GF != NULL);
  return FpElement(_GF, _GF->sub(this->_val, other._val));
}

template<typename R, typename S>
FpElement<R,S> FpElement<R,S>::operator*(const FpElement<R, S> &other) const
{
  assert(_GF != NULL);
  return FpElement(_GF, _GF->mul(this->_val, other._val));
}

template<typename R, typename S>
FpElement<R,S> FpElement<R,S>::operator/(const FpElement<R, S> &other) const
{
  assert((_GF != NULL) && (other != 0));
  R elt = _GF->mod(other._val).lift();
  return FpElement(_GF, _GF->mul(this->_val, _GF->inverse(elt)));
}

template<typename R, typename S>
FpElement<R,S> FpElement<R,S>::inverse(void) const
{
  assert((_GF != NULL) && ((*this) != 0));
  R elt = _GF->mod(this->_val).lift();
  return FpElement(_GF, _GF->inverse(elt));
}

template<typename R, typename S>
FpElement<R,S> & FpElement<R,S>::operator+=(const FpElement<R, S> &other)
{
  assert(_GF != NULL);
  _val = _GF->add(this->_val, other._val); return (*this);
}

template<typename R, typename S>
FpElement<R,S> & FpElement<R,S>::operator-=(const FpElement<R, S> &other)
{
  assert(_GF != NULL);
  _val = _GF->sub(this->_val, other._val); return (*this);
}

template<typename R, typename S>
FpElement<R,S> & FpElement<R,S>::operator*=(const FpElement<R, S> &other)
{
  assert(_GF != NULL);
  _val = _GF->mul(this->_val, other._val); return (*this);
}

template<typename R, typename S>
FpElement<R,S> & FpElement<R,S>::operator/=(const FpElement<R, S> &other)
{
  assert((_GF != NULL) && (other != 0));
  R elt = _GF->mod(other._val).lift();
  _val = _GF->mul(this->_val, _GF->inverse(elt)); return (*this);
}

template<typename R, typename S>
FpElement<R,S> FpElement<R,S>::sqrt() const
{
  assert(_GF != NULL);
  return FpElement(_GF, _GF->sqrt(this->_val));
}

// assignment and conversion
template<typename R, typename S>
FpElement<R,S> & FpElement<R,S>::operator=(const FpElement<R, S> &other) 
{
  if (this != (&other)) {
    this->_GF = other._GF;
    this->_val = other._val;
  }
  return (*this);
}

  
//boolean
template<typename R, typename S>
bool FpElement<R,S>::operator==(const FpElement<R, S> &other) const {
  assert( (_GF != 0) && (other._GF != 0) ); 
  if (this->_GF->prime() != other._GF->prime()) return false;
  return ((this->_val - other._val) % (this->_GF->prime()) == 0);
}

template<typename R, typename S>
bool FpElement<R,S>::operator==(const R &other) const {
  assert(_GF != 0);
  return ((this->_val - other) % (this->_GF->prime()) == 0);
}

template<typename R, typename S>
bool FpElement<R,S>::operator!=(const R &other) const {
  assert(_GF != 0);
  return ((this->_val - other) % (this->_GF->prime()) != 0);
}

template<typename R, typename S>
bool FpElement<R,S>::isSquare(void) const {
  assert(_GF != 0);
  return ((this->_GF->legendre(this->_val)) >= 0);
}
