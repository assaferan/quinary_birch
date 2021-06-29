#include "COMB_Stor.h"

void VSLT_ZeroAssignment(IntegerSubsetStorage & VSLT)
{
  int MaxElement=VSLT.MaxElement;
  int Maxp2=MaxElement+2;
  for (int iVert=0; iVert<Maxp2; iVert++) {
    VSLT.ListNext[iVert]=-1;
    VSLT.ListPrev[iVert]=-1;
  }
  VSLT.ListNext[MaxElement  ] = MaxElement+1;
  VSLT.ListPrev[MaxElement+1] = MaxElement;
}


IntegerSubsetStorage VSLT_InitializeStorage(int const& MaxElement)
{
  IntegerSubsetStorage VSLT;
  int Maxp2 = MaxElement+2;
  std::vector<int> ListNext(Maxp2);
  std::vector<int> ListPrev(Maxp2);
  VSLT.ListNext = ListNext;
  VSLT.ListPrev = ListPrev;
  VSLT.MaxElement = MaxElement;
  VSLT_ZeroAssignment(VSLT);
  return VSLT;
}


int VSLT_NrElement(IntegerSubsetStorage const& VSLT)
{
  int MaxElt=VSLT.MaxElement+1;
  int pos=MaxElt-1;
  int NbElt=0;
  while(true) {
    int posNext=VSLT.ListNext[pos];
    if (posNext == MaxElt)
      return NbElt;
    pos=posNext;
    NbElt++;
  }
}

int VSLT_TheFirstPosition(IntegerSubsetStorage const& VSLT)
{
  return VSLT.ListNext[VSLT.MaxElement];
}

bool VSLT_IsItInSubset(IntegerSubsetStorage const& VSLT, int const& pos)
{
  return VSLT.ListNext[pos] != -1;
}

void VSLT_StoreValue(IntegerSubsetStorage & VSLT, int const& pos)
{
  int posAfter=VSLT.ListNext[VSLT.MaxElement];
  VSLT.ListNext[VSLT.MaxElement] = pos;
  VSLT.ListNext[pos] = posAfter;
  VSLT.ListPrev[posAfter] = pos;
  VSLT.ListPrev[pos] = VSLT.MaxElement;
}



void VSLT_RemoveValue(IntegerSubsetStorage & VSLT, int const& pos)
{
  int posNext=VSLT.ListNext[pos];
  int posPrev=VSLT.ListPrev[pos];
  VSLT.ListNext[posPrev] = posNext;
  VSLT.ListPrev[posNext] = posPrev;
  VSLT.ListNext[pos] = -1;
  VSLT.ListPrev[pos] = -1;
}


bool VSLT_IsEmpty(IntegerSubsetStorage const& VSLT)
{
  return VSLT.ListNext[VSLT.MaxElement] == VSLT.MaxElement+1;
}

EncapsDynamicBitset::EncapsDynamicBitset(std::size_t const& eMax) : eVect(eMax)
{
  TotalCount=0;
}

EncapsDynamicBitset::~EncapsDynamicBitset()
{
}

std::size_t EncapsDynamicBitset::count() const
{
  return TotalCount;
}

std::size_t EncapsDynamicBitset::find_first() const
{
  return eVect.find_first();
}

void EncapsDynamicBitset::set(std::size_t const& pos, bool const& eVal)
{
  if (eVect[pos] && !eVal)
    TotalCount--;
  if (!eVect[pos] && eVal)
    TotalCount++;
  eVect[pos]=eVal;
}

bool EncapsDynamicBitset::get(std::size_t const& pos) const
{
  return eVect[pos];
}

bool EncapsDynamicBitset::empty() const
{
  if (TotalCount == 0)
    return true;
  return false;
}
