#include "TypeConversion.h"

int IntFloor(double const& x)
{
  return int(floor(x));
}

void NearestInteger_double_int(double const& xI, int & xO)
{
  //  std::cerr << "Temp_common : NearestInteger\n";
  double xRnd_d=round(xI);
  int xRnd_z=int(xRnd_d);
  //  std::cerr << "xI=" << xI << "\n";
  auto GetErr=[&](int const& u) -> double {
    double diff = double(u) - xI;
    return T_abs(diff);
  };
  double err=GetErr(xRnd_z);
  //  std::cerr << "err=" << err << "\n";
  while(true) {
    bool IsOK=true;
    for (int i=0; i<2; i++) {
      int shift=2*i -1;
      int xTest = xRnd_z + shift;
      double TheErr=GetErr(xTest);
      //      std::cerr << "i=" << i << " shift=" << shift << " xTest=" << xTest << " TheErr=" << TheErr << "\n";
      if (TheErr < err) {
	IsOK=false;
	xRnd_z=xTest;
      }
    }
    if (IsOK)
      break;
  }
  xO=xRnd_z;
}

void TYPE_CONVERSION(double const& a1, double & a2)
{
  a2 = a1;
}

void TYPE_CONVERSION(double const& a1, uint8_t & a2)
{
  a2 = uint8_t(a1);
}

void TYPE_CONVERSION(double const& a1, int & a2)
{
  a2 = int(a1);
}

void TYPE_CONVERSION(double const& a1, long & a2)
{
  a2 = long(a1);
}

void TYPE_CONVERSION(int const& a1, double & a2)
{
  a2 = double(a1);
}

void TYPE_CONVERSION(uint8_t const& a1, double & a2)
{
  a2 = double(a1);
}

void TYPE_CONVERSION(long const& a1, double & a2)
{
  a2 = double(a1);
}



