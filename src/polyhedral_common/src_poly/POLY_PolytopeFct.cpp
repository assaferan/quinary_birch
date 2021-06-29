#include "POLY_PolytopeFct.h"

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
