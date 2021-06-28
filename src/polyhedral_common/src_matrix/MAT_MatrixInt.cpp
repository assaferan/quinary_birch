#include "MAT_MatrixInt.cpp"

int IsVectorPrimitive(MyVector<int> const& TheV)
{
  int n=TheV.size();
  int TheGCD=TheV(0);
  for (int i=1; i<n; i++) {
    int eValI=TheV(i);
    GCD_int<int> eRec=ComputePairGcd(TheGCD, eValI);
    TheGCD=eRec.gcd;
  }
  if (abs(TheGCD) == 1)
    return 1;
  return 0;
}
