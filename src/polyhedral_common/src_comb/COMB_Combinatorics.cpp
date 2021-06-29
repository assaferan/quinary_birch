#include "COMB_Combinatorics.h"

int BlockIteration::IncrementShow()
{
  for (int i=0; i<dim; i++)
    if (eVect[i] < size-1) {
      eVect[i]++;
      for (int j=0; j<i; j++)
	eVect[j]=0;
      return i;
    }
  return -1;
}

void BlockIteration::IncrementSilent()
{
  for (int i=0; i<dim; i++)
    if (eVect[i] < size-1) {
      eVect[i]++;
      for (int j=0; j<i; j++)
	eVect[j]=0;
      return;
    }
  std::cerr << "Should not reach that stage\n";
  throw TerminalException{1};
}

std::vector<int> BlockIteration::GetVect() const
{
  return eVect;
}

int BlockIteration::GetNbPoss() const
{
  int eRet=1;
  for (int i=0; i<dim; i++)
    eRet *= size;
  return eRet;
}

void BlockCppIterator::IteratorContain::single_increase()
{
  for (int i=0; i<dim_iter; i++)
    if (U[i] < size_iter-1) {
      U[i]++;
      for (int j=0; j<i; j++)
	U[j]=0;
      return;
    }
  U = {};
}

std::vector<int> const& BlockCppIterator::IteratorContain::operator*()
{
  return U;
}

BlockCppIterator::IteratorContain & BlockCppIterator::IteratorContain::operator++()
{
  single_increase();
  return *this;
}

BlockCppIterator::IteratorContain BlockCppIterator::IteratorContain::operator++(int)
{
  BlockCppIterator::IteratorContain tmp = *this;
  single_increase();
  return tmp;
}

bool BlockCppIterator::IteratorContain::operator!=(IteratorContain const& iter)
{
  if (iter.dim_iter != dim_iter)
    return true;
  if (iter.size_iter != size_iter)
    return true;
  if (iter.U.size() != U.size())
    return true;
  for (size_t i=0; i<U.size(); i++)
    if (iter.U[i] != U[1])
      return true;
  return false;
}

bool BlockCppIterator::IteratorContain::operator==(IteratorContain const& iter)
{
  if (iter.dim_iter != dim_iter)
    return false;
  if (iter.size_iter != size_iter)
    return false;
  if (iter.U.size() != U.size())
    return false;
  for (size_t i=0; i<U.size(); i++)
    if (iter.U[i] != U[1])
      return false;
  return true;
}

// The iterator business

BlockCppIterator::IteratorContain BlockCppIterator::cbegin() const
{
  return { dim, size, std::vector<int>(dim,0) };
}

BlockCppIterator::IteratorContain BlockCppIterator::cend() const
{
  return {dim, size, {}};
}

BlockCppIterator::IteratorContain BlockCppIterator::begin() const
{
  return { dim, size, std::vector<int>(dim,0) };
}

BlockCppIterator::IteratorContain BlockCppIterator::end() const
{
  return {dim, size, {}};
}

int BlockIterationMultiple::IncrementShow()
{
  for (size_t i=0; i<dim; i++)
    if (eVect[i] < ListSize[i]-1) {
      eVect[i]++;
      for (size_t j=0; j<i; j++)
	eVect[j]=0;
      return i;
    }
  return -1;
}

void BlockIterationMultiple::IncrementSilent()
{
  for (size_t i=0; i<dim; i++)
    if (eVect[i] < ListSize[i]-1) {
      eVect[i]++;
      for (size_t j=0; j<i; j++)
	eVect[j]=0;
      return;
    }
  std::cerr << "Should not reach that stage\n";
  throw TerminalException{1};
}

std::vector<int> BlockIterationMultiple::GetVect() const
{
  return eVect;
}

size_t BlockIterationMultiple::GetNbPoss() const
{
  size_t eRet=1;
  for (size_t i=0; i<dim; i++)
    eRet *= ListSize[i];
  return eRet;
}

std::vector<int> RandomPermutation(int const& n)
{
  std::vector<int> RetList(n,-1);
  for (int i=0; i<n; i++) {
    int rnd_pos = rand() % (n - i);
    int idx=0;
    bool IsAssigned = false;
    for (int j=0; j<n; j++) {
      if (RetList[j] == -1 && !IsAssigned) {
        if (rnd_pos == idx) {
          RetList[j] = i;
          IsAssigned = true;
        }
        idx++;
      }
    }
  }
  return RetList;
}

int BlockIteration01::IncrementShow()
{
  for (int i=0; i<dim; i++)
    if (eFace[i] == 0) {
      eFace[i]=1;
      for (int j=0; j<i; j++)
	eFace[j]=0;
      return i;
    }
  return -1;
}

void BlockIteration01::IncrementSilent()
{
  for (int i=0; i<dim; i++)
    if (eFace[i] == 0) {
      eFace[i]=1;
      for (int j=0; j<i; j++)
	eFace[j]=0;
      return;
    }
  std::cerr << "Should not reach that stage\n";
  throw TerminalException{1};
}

Face BlockIteration01::GetFace() const
{
  return eFace;
}

int BlockIteration01::GetNbPoss() const
{
  int eRet=1;
  for (int i=0; i<dim; i++)
    eRet *= 2;
  return eRet;
}
