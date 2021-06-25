// std includes

#include <sstream>

// flint/antic includes

#include "antic/nf.h"
#include "antic/nf_elem.h"

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

  nf_t nf;
  nf_elem_t a,b,c;
  fmpq_t atrace, btrace, ctrace, ctrace2;

  flint_rand_t state;

  flint_randinit(state);
  
  nf_elem_init(a, nf);
  nf_elem_init(b, nf);
  nf_elem_init(c, nf);

  fmpq_init(atrace);
  fmpq_init(btrace);
  fmpq_init(ctrace);
  fmpq_init(ctrace2);

  nf_elem_add(c, a, b, nf);
  nf_elem_trace(atrace, a, nf);
  nf_elem_trace(btrace, b, nf);
  nf_elem_trace(ctrace, c, nf);
  fmpq_add(ctrace2, atrace, btrace);

  bool result = (fmpq_equal(ctrace, ctrace2));

  fmpq_poly_print_pretty(nf->pol, "x"); std::cerr << std::endl;
  nf_elem_print_pretty(a, nf, "x"); std::cerr << std::endl;
  nf_elem_print_pretty(b, nf, "x"); std::cerr << std::endl;
  nf_elem_print_pretty(c, nf, "x"); std::cerr << std::endl;
  fmpq_print(atrace); std::cerr << std::endl;
  fmpq_print(btrace); std::cerr << std::endl;
  fmpq_print(ctrace); std::cerr << std::endl;
  fmpq_print(ctrace2); std::cerr << std::endl;

  fmpq_clear(atrace);
  fmpq_clear(btrace);
  fmpq_clear(ctrace);
  fmpq_clear(ctrace2);

  nf_elem_clear(a, nf);
  nf_elem_clear(b, nf);
  nf_elem_clear(c, nf);

  nf_clear(nf);

  flint_randclear(state);
  flint_cleanup();

  runBirchTests(num_evs);
  
  return 0;
}
