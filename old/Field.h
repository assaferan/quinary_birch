#ifndef __FIELD_H
#define __FIELD_H

#include "Ring.h"

// Figure out what should be here, if any
class Field : public Ring
{
public:
};

class FieldElement : public RingElement
{
 public:
  
  // arithmetic
  virtual FieldElement operator/(const FieldElement & ) const = 0;

  virtual FieldElement & operator/=(const FieldElement & ) = 0;

  virtual FieldElement inverse(void) const = 0;

};

#endif // __FIELD_H
