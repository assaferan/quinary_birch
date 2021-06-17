#ifndef __EIGENVECTOR_H_
#define __EIGENVECTOR_H_

#include <algorithm>
#include "NumberFieldElement.h"
#include "SetCover.h"

template<typename R>
class Eigenvector
{
public:
  Eigenvector() = default;
  Eigenvector(const Eigenvector<R>& other) = default;
  Eigenvector(Eigenvector<R>&& other) = default;

  Eigenvector<R>& operator=(const Eigenvector<R>& other) = default;
  Eigenvector<R>& operator=(Eigenvector<R>&& other) = default;

  Eigenvector(std::vector< NumberFieldElement<Z> >&& data, W64 conductor_index)
  {
    this->_data = std::vector< NumberFieldElement<Z> >(data);
    this->_conductor_index = conductor_index;
  }

  inline const std::vector< NumberFieldElement<Z> >& data(void) const
  { return this->_data; }

  inline size_t size(void) const
  { return this->_data.size(); }

  inline W64 conductorIndex(void) const
  { return this->_conductor_index; }

  inline void repIndex(size_t pos)
  { this->_rep_index = pos; }

  inline size_t repIndex(void) const
  { return this->_rep_index; }

  inline NumberFieldElement<Z> operator[](size_t pos) const
  { return this->_data[pos]; }

private:
  std::vector< NumberFieldElement<Z> > _data;
  W64 _conductor_index;
  size_t _rep_index;
};

template<typename R, size_t n>
class EigenvectorManager
{
  friend class Genus<R,n>;

public:
  EigenvectorManager() = default;

  void addEigenvector(const Eigenvector<R>& vector);

  inline size_t size(void) const
  { return this->_eigenvectors.size(); }

  void finalize(void);

  inline const Eigenvector<R>& operator[](size_t index) const
  { return this->_eigenvectors[index]; }

private:
  bool _finalized = false;
  size_t _dimension = 0;
  size_t _stride = 0;
  std::vector<Eigenvector<R>> _eigenvectors;
  std::vector<Z32> _strided_eigenvectors;
  std::vector<W64> _conductors;
  std::vector<std::vector<Z64>> _position_lut;
  std::vector<Z64> _indices;
  W64 _conductor_primes;
};

#include "Eigenvector.inl"

#endif // __EIGENVECTOR_H_
