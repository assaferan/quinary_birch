#include "POLY_PolytopeFct.h"

std::vector<int> Dynamic_bitset_to_vectorint(Face const& eList)
{
  int nb=eList.count();
  int aRow=eList.find_first();
  std::vector<int> retList(nb);
  for (int i=0; i<nb; i++) {
    retList[i] = aRow;
    aRow=eList.find_next(aRow);
  }
  return retList;
}

void PrintListOrbit(std::ostream &os, std::vector<Face> const& ListOrbit)
{
  int iOrbit, nbOrbit, siz, i, eVal;
  nbOrbit=ListOrbit.size();
  os << "nbOrbit=" << nbOrbit << "\n";
  for (iOrbit=0; iOrbit<nbOrbit; iOrbit++) {
    Face eInc=ListOrbit[iOrbit];
    siz=eInc.count();
    os << "O" << iOrbit+1 << ": inc=" << siz << "list=";
    eVal=eInc.find_first();
    for (i=0; i<siz; i++) {
      os << " " << eVal;
      eVal=eInc.find_next(eVal);
    }
    os << "\n";
  }
}
