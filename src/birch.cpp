// std includes

#include <sstream>

#include "polyhedral_common/src_latt/MatrixCanonicalForm.h"
#include "TestBirch.h"


int main(int argc, char* argv[])
{
  size_t num_evs = 0;
  
  if (argc > 2) {
    std::cerr << "Usage: " << argv[0] << " [num_evs]" << std::endl;
    std::cerr << "where [num_evs] is a number specifying how many eigenvalues to test." << std::endl;
  }

  if (argc == 2) {
    std::stringstream num_evs_str(argv[1]);
    num_evs_str >> num_evs;
  }

  std::ifstream is("src/polyhedral_common/Examples/ConwaySloane_LatticeNoBasis");
  MyMatrix<Z32> eMat = ReadMatrix<Z32>(is);
  Canonic_PosDef<Z32,Z32> RetF = ComputeCanonicalForm<Z32,Z32>(eMat);
  
  runBirchTests(num_evs);
  
  return 0;
}
