#include <functional>
#include <iostream>

#include "polyhedral_common/src_basic/ExceptionEnding.h"

#include "factorizations.h"
#include "NumberTheory.h"

std::vector<mpz_class> FactorsInt(mpz_class const& x)
{
  if (!IsInteger(x)) {
    std::cerr << "FactorsInt requires the input to be integral\n";
    throw TerminalException{1};
  }
  std::vector<mpz_class> ListDiv;
  mpz_class TheWork=x;
  std::function<void(void)> GetOneDivisor=[&]() -> void {
    mpz_class eVal=2;
    while(true) {
      mpz_class TheRes=TheWork % eVal;
      if (TheRes == 0) {
	ListDiv.push_back(eVal);
	TheWork = TheWork / eVal;
	return;
      }
      eVal += 1;
    }
  };
  while(true) {
    if (TheWork == 1)
      break;
    GetOneDivisor();
  }
  return ListDiv;
}

std::vector<mpq_class> FactorsInt(mpq_class const& x)
{
  mpz_class x_z=x.get_num();
  std::vector<mpz_class> LFact=FactorsInt(x_z);
  std::vector<mpq_class> LFact_q;
  for (auto & eVal : LFact) {
    mpq_class eVal_q=eVal;
    LFact_q.push_back(eVal_q);
  }
  return LFact_q;
}
