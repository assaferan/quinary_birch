// std includes

#include <sstream>

#include "TestBirch.h"

bool handleFlagInt(const std::string & flag_name, const std::string & param_str, int & flag_val);
int printParamDesc(char* argv[]);

int main(int argc, char* argv[])
{
  bool do_tests = false;
  int num_evs = 0;
  ReductionMethod alg = GREEDY_FULL;
  int disc = 0;
  
  if (argc > 5) {
    return printParamDesc(argv);
  }

  for (int i = 1; i < argc; i++) {
    std::string param_str(argv[i]);
    bool is_valid = false;
    bool is_disc, is_nevs;
    if (param_str == "-tests") {
      do_tests = true;
      is_valid = true;
    }

    is_disc = handleFlagInt("d", param_str, disc);
    is_valid = is_valid || is_disc;
    is_nevs = handleFlagInt("nevs", param_str, num_evs);
    is_valid = is_valid || is_nevs;
    
    if (param_str.substr(0,5) == "-red=") {
      param_str = param_str.substr(5, param_str.length()-5);
      if (param_str == "canonical") {
	alg = CANONICAL_FORM;
	is_valid = true;
      }
      else if (param_str == "greedy") {
	alg = GREEDY;
	is_valid = true;
      }
      if (param_str == "greedy_full") {
	alg = GREEDY_FULL;
	is_valid = true;
      }
      if (param_str == "minkowski") {
	alg = MINKOWSKI;
	is_valid = true;
      }
    }
    
    if (!is_valid)
      return printParamDesc(argv);
  }

  if (disc > 0)
    runQuinaryBirch(disc,num_evs,alg);
  if (do_tests)
    runBirchTests(num_evs,alg);
  
  return 0;
}

bool handleFlagInt(const std::string & flag_name, const std::string & param_str, int & flag_val)
{
  size_t param_len = param_str.length();

  std::string full_flag_name = "-" + flag_name + "=";

  size_t flag_len = full_flag_name.length();
  
  if (param_str.substr(0,flag_len) == full_flag_name) {
    std::stringstream flag_str(param_str.substr(flag_len,param_len-flag_len));
    flag_str >> flag_val;
    return true;
  }

  return false;
}

int printParamDesc(char* argv[])
{
  std::cerr << "Usage: " << argv[0] << " [-tests] [-d=disc] [-nevs=num_evs] [-red=reduction_method]" << std::endl;
  std::cerr << "[disc] is the discriminant to compute eigenvalues for," << std::endl;
  std::cerr << "[num_evs] is a number specifying how many eigenvalues to test," << std::endl;
  std::cerr << "[reduction_method] is either \"greedy\", \"greedy_full\" or \"canonical\", the default being greedy." << std::endl;
  std::cerr << "If the flag -tests is supplied, runs standard tests, with the parameters num_evs, reduction_method as supplied." << std::endl;
  return -1;
}
