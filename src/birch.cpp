// std includes

#include <sstream>

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

  runBirchTests(num_evs);
  
  return 0;
}
