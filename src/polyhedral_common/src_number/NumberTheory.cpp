#include "NumberTheory.h"

mpz_class GetRingElement(mpq_class const& eVal)
{
  return eVal.get_num();
}

namespace std {
  std::string to_string(const mpz_class& e_val)
  {
    std::stringstream s;
    s << e_val;
    std::string converted(s.str());
    return converted;
  };
  std::string to_string(const mpq_class& e_val)
  {
    std::stringstream s;
    s << e_val;
    std::string converted(s.str());
    return converted;
  };
}
