// std includes

#include <sstream>

#include "TestBirch.h"

int printParamDesc(char* argv[]);

int main(int argc, char* argv[])
{
  size_t num_evs = 0;
  ReductionMethod alg = GREEDY;
  Z64 disc;
  
  if ((argc > 4) || (argc < 2) ){
    return printParamDesc(argv);
  }

  std::stringstream disc_str(argv[1]);
  disc_str >> disc;

  if (argc >= 3) {
    std::stringstream num_evs_str(argv[2]);
    num_evs_str >> num_evs;
  }

  if (argc == 4) {
    std::string alg_str(argv[3]);
    int cmp_canonical = alg_str.compare("canonical");
    if (cmp_canonical == 0) {
      alg = CANONICAL_FORM;
    }
    int cmp_greedy = alg_str.compare("greedy");
    if (cmp_greedy == 0) {
      alg = GREEDY;
    }
    if ((cmp_canonical != 0) && (cmp_greedy != 0)) {
      return printParamDesc(argv);
    }
  }
  
  runQuinaryBirch(disc,num_evs,alg);
  runBirchTests(num_evs,alg);
  
  return 0;
}

int printParamDesc(char* argv[])
{
  std::cerr << "Usage: " << argv[0] << " [disc] [num_evs] [reduction_method]" << std::endl;
  std::cerr << "where [disc] is the discriminant to compute eigenvalues for," << std::endl;
  std::cerr << "[num_evs] is a number specifying how many eigenvalues to test," << std::endl;
  std::cerr << "[reduction_method] is either \"greedy\" or \"canonical\", the default being greedy." << std::endl;
  return -1;
}
