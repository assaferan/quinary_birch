#include <cstdlib>
#include "ExceptionEnding.h"

void TerminalEnding()
{
#ifdef DEBUG
  int val1=4;
#endif
  assert(val1 == 5);
  std::exit(1);
}
