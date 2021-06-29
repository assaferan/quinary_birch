#ifndef INCLUDE_COMBINATORICS_ELEM
#define INCLUDE_COMBINATORICS_ELEM

#include "polyhedral_common/src_basic/Temp_common.h"
#include "polyhedral_common/src_matrix/MAT_Matrix.h"

std::vector<int> BinomialStdvect_First(int const& k);
bool BinomialStdvect_Increment(int const&n, int const&k, std::vector<int> & Tvect);
MyMatrix<int> BuildSet(int const& n, int const& Nval);
int PositionBuildSet(int const& n, int const& Nval, MyVector<int> const& V);

#endif


