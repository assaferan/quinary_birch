#ifndef INCLUDE_FACE_BITSET
#define INCLUDE_FACE_BITSET

// Boost libraries

#include "polyhedral_common/src_basic/Temp_common.h"
#include <bitset>
#include <boost/dynamic_bitset.hpp>

typedef boost::dynamic_bitset<> Face;
std::vector<int> FaceToVector(Face const& eSet);
std::vector<int> FaceTo01vector(Face const& eSet);
void WriteFace(std::ostream & os, Face const& eList);
Face ReadFace(std::istream & is);
std::vector<Face> ReadListFace(std::istream & is);
void WriteListFace(std::ostream & os, std::vector<Face> const& ListFace);
void WriteFaceGAP(std::ostream &os, Face const& f);
void WriteListFaceGAP(std::ostream & os, std::vector<Face> const& ListFace);
void WriteListFaceGAPfile(std::string const& eFile, std::vector<Face> const& ListFace);
// We require x and y to be of the same size
bool operator<(Face const& x, Face const& y);
void PrintVectInt(std::ostream &os, Face const& eList);
Face FullFace(int const& len);
ulong FaceToUnsignedLong(Face const& f);
Face UnsignedLongToFace(int const& len, ulong const& eVal);
void VectVectInt_Magma_Print(std::ostream &os, std::vector<Face> const&ListOrbit);

#endif
