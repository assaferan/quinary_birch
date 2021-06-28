#include "NumberTheory.h"

mpz_class GetRingElement(mpq_class const& eVal)
{
  return eVal.get_num();
}
