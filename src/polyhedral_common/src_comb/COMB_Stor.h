#ifndef POLYTOPE_STOR
#define POLYTOPE_STOR

#include "polyhedral_common/src_basic/Temp_common.h"
#include "Boost_bitset.h"

//
// This is a set of functionality for storing bits
// with different complexity theories and memory requirements.
// Those classes are supposed to be used as template arguments.
//
// Each class implements following subset of dynamic_subset:
// a[i] = val with val a boolean, i.e. 0/1 or true/false.
// val  = a[i]
// find_first()  : first true element of the list
// empty()       : true if empty, false otherwise
// count()       : number of true entries
//
// Possible available classes:
// dynamic_bitset (aliases to Face)
//  ---find_first  : O(n)
//  ---data access : O(1)
//  ---count       : O(n)
// Memory expenses : 2^n + overhead
//
// DoubleList:
//  ---find_first  : O(1)
//  ---data access : O(1)
//  ---count       : O(n)
// Memory expenses : 2n * 2^n + overhead
//

struct IntegerSubsetStorage {
  int MaxElement;
  std::vector<int> ListNext;
  std::vector<int> ListPrev;
};


void VSLT_ZeroAssignment(IntegerSubsetStorage & VSLT);
IntegerSubsetStorage VSLT_InitializeStorage(int const& MaxElement);
int VSLT_NrElement(IntegerSubsetStorage const& VSLT);
int VSLT_TheFirstPosition(IntegerSubsetStorage const& VSLT);
bool VSLT_IsItInSubset(IntegerSubsetStorage const& VSLT, int const& pos);
void VSLT_StoreValue(IntegerSubsetStorage & VSLT, int const& pos);
void VSLT_RemoveValue(IntegerSubsetStorage & VSLT, int const& pos);
bool VSLT_IsEmpty(IntegerSubsetStorage const& VSLT);

// This is an artificial class for allowing the operation
// a[i]=a
template<typename Tint>
struct DoubleList {
public:
  DoubleList(const DoubleList&) = delete;
  DoubleList& operator=(const DoubleList&) = delete;
  DoubleList(DoubleList&&) = delete;
  DoubleList() = delete;
  DoubleList(Tint const& eMax) : MaxElement(eMax)
  {
    Tint Maxp2=eMax + 2;
    ListNext = std::vector<Tint>(Maxp2, Maxp2);
    ListPrev = std::vector<Tint>(Maxp2, Maxp2);
    ListNext[MaxElement]=MaxElement+1;
    ListPrev[MaxElement+1]=MaxElement;
  }
  ~DoubleList()
  {
  }
  Tint count() const
  {
    Tint Maxp1=MaxElement+1;
    Tint pos=MaxElement;
    Tint NbElt=0;
    //    std::cerr << "MaxElement=" << MaxElement << " Maxp1=" << Maxp1 << "\n";
    while(true) {
      //      std::cerr <<  "   pos=" << pos << " ListNext[pos]=" << ListNext[pos] << "\n";
      pos=ListNext[pos];
      if (pos == Maxp1)
	return NbElt;
      NbElt++;
    }
    std::cerr << "We should not reach that stage\n";
    throw TerminalException{1};
  }
  bool empty() const
  {
    Tint Maxp1=MaxElement+1;
    //    std::cerr << "MaxElement=" << MaxElement << " Maxp1=" << Maxp1 << " ListNext[]=" << ListNext[MaxElement] << "\n";
    if (ListNext[MaxElement] == Maxp1)
      return true;
    return false;
  }
  Tint find_first() const
  {
    return ListNext[MaxElement];
  }
  void set(Tint const& pos, bool const& eVal)
  {
    Tint Maxp2=MaxElement+2;
    if (eVal) {
      if (ListNext[pos] == Maxp2) {
	Tint posAfter=ListNext[MaxElement];
	ListNext[MaxElement]=pos;
	ListNext[pos]=posAfter;
	ListPrev[posAfter]=pos;
	ListPrev[pos]=MaxElement;
      }
    }
    else {
      if (ListNext[pos] != Maxp2) {
	Tint posNext=ListNext[pos];
	Tint posPrev=ListPrev[pos];
	ListNext[posPrev]=posNext;
	ListPrev[posNext]=posPrev;
	ListNext[pos]=Maxp2;
	ListPrev[pos]=Maxp2;
      }
    }
  }
  bool get(Tint const& pos) const
  {
    Tint Maxp2=MaxElement+2;
    if (ListNext[pos] == Maxp2)
      return false;
    return true;
  }
private:
  Tint MaxElement;
  std::vector<Tint> ListNext;
  std::vector<Tint> ListPrev;
};


struct EncapsDynamicBitset {
public:
  EncapsDynamicBitset(const EncapsDynamicBitset&) = delete;
  EncapsDynamicBitset& operator=(const EncapsDynamicBitset&) = delete;
  EncapsDynamicBitset(EncapsDynamicBitset&&) = delete;
  EncapsDynamicBitset() = delete;
  EncapsDynamicBitset(std::size_t const& eMax);
  ~EncapsDynamicBitset();
  
  std::size_t count() const;
  std::size_t find_first() const;
  void set(std::size_t const& pos, bool const& eVal);
  bool get(std::size_t const& pos) const;
  bool empty() const;
  
private:
  std::size_t TotalCount;
  boost::dynamic_bitset<> eVect;
};



template<typename T>
struct StackStorage {
public:
  StackStorage()
  {
    TheSize=0;
  }
  void push_back(T const& val)
  {
    if (TheSize == ListElt.size()) {
      TheSize++;
      ListElt.push_back(val);
    }
    else {
      ListElt[TheSize]=val;
      TheSize++;
    }
  }
  T pop()
  {
    T val=ListElt[TheSize-1];
    TheSize--;
    return val;
  }
  size_t size()
  {
    return TheSize;
  }
private:
  std::vector<T> ListElt;
  size_t TheSize;
};





#endif
