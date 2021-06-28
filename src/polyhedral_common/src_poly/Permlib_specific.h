#ifndef INCLUDE_PERMLIB_SPECIFIC_H
#define INCLUDE_PERMLIB_SPECIFIC_H


#include "polyhedral_common/src_comb/Boost_bitset.h"
#include "polyhedral_common/src_basic/hash_functions.h"
#include <permlib/permlib_api.h>


int OnPoints(int const& i, permlib::Permutation const& elt);

//
// Specific implementation with permlib
//

typedef std::shared_ptr<permlib::PermutationGroup> PermutationGroupPtr;
typedef boost::dynamic_bitset<> DsetList;

/*
std::vector<permlib::dom_int> GetBaseGroup(TheGroupFormat const& eGRP);
*/
std::set<int> GetSetFrom_DB(Face const& eList);
Face Face_EquivSimplification(Face const& eFace);
std::pair<bool,permlib::Permutation> PERMLIB_TestEquivalenceGeneralNaked(int const& n, PermutationGroupPtr const& group, Face const& eList1, Face const& eList2, int const& eMethod);

template<typename Tint>
std::pair<bool,permlib::Permutation> PERMLIB_TestEquivalenceGeneral(int const& n, PermutationGroupPtr const& group, Tint const& grp_size, Face const& eList1, Face const& eList2)
{
  Tint MaxSize=10000;
  int eMethod = 1;
  if (grp_size < MaxSize)
    eMethod=0;
  return PERMLIB_TestEquivalenceGeneralNaked(n, group, eList1, eList2, eMethod);
}

template<typename Tint>
bool PERMLIB_TestEquivalence(int const& n, PermutationGroupPtr const& group, Tint const& grp_size, Face const& eList1, Face const& eList2)
{
  return PERMLIB_TestEquivalenceGeneral(n, group, grp_size, eList1, eList2).first;
}

Face PERMLIB_Canonicalization(int const& n, PermutationGroupPtr const& group, Face const& eList);
PermutationGroupPtr PERMLIB_GetStabilizer_general(PermutationGroupPtr const& group, Face const& eList, int const& opt);

template<typename Tint>
PermutationGroupPtr PERMLIB_GetStabilizer(PermutationGroupPtr const& group, Tint const& grp_size, Face const& eList)
{
  Tint MaxSize=10000;
  if (grp_size < MaxSize) {
    return PERMLIB_GetStabilizer_general(group, eList, 0);
  }
  else {
    return PERMLIB_GetStabilizer_general(group, eList, 1);
  }
}

template<typename Tint_inp>
struct TheGroupFormat {
private:
  int n;
  Tint_inp e_size;
  PermutationGroupPtr group;
  TheGroupFormat(int const& _n, Tint_inp const& _size, PermutationGroupPtr const& _group) : n(_n), e_size(_size), group(_group)
  {
  }
public:
  using Tint = Tint_inp;
  using Telt = permlib::Permutation;
  TheGroupFormat(std::vector<permlib::Permutation> const& ListPerm, int const& n_inp)
  {
    std::vector<permlib::Permutation::ptr> generatorList;
    for (auto & eGen : ListPerm) {
      std::vector<int> v(n);
      for (int i=0; i<n; i++)
        v[i]=eGen.at(i);
      generatorList.push_back(permlib::Permutation::ptr(new permlib::Permutation(v)));
    }
    n = n_inp;
    group = construct(n, generatorList.begin(), generatorList.end());
    e_size = group->order<Tint>();
  }
  TheGroupFormat(int const& n_inp) : TheGroupFormat({}, n) 
  {
  }
  TheGroupFormat() : TheGroupFormat(0)
  {
  }
  Face CanonicalImage(Face const& eFace) const
  {
    return PERMLIB_Canonicalization(n, group, eFace);
  }
  TheGroupFormat<Tint> Stabilizer_OnSets(Face const& f) const
  {
    PermutationGroupPtr group_stab = PERMLIB_GetStabilizer(group, e_size, f);
    Tint_inp size_stab = group_stab->order<Tint_inp>();
    return TheGroupFormat(n, size_stab, group_stab);
  }
  std::pair<bool,permlib::Permutation> RepresentativeAction_OnSets(Face const& f1, Face const& f2) const
  {
    return PERMLIB_TestEquivalenceGeneral(n, group, e_size, f1, f2);
  }
  std::vector<Telt> GeneratorsOfGroup() const
  {
    // copy operation, but that is what it is.
    std::vector<Telt> LGen;
    for (auto & eval : group->S) {
      LGen.push_back(*eval);
    }
    return LGen;
  }
  Tint size() const
  {
    return e_size;
  }
  int n_act() const
  {
    return n;
  }
  PermutationGroupPtr get_group() const // USe only by some function specific to permlib.
  {
    return group;
  }
};

namespace std {
  template<>
  struct hash<permlib::Permutation>
  {
    std::size_t operator()(const permlib::Permutation & eElt) const;
  };
}

void WriteVectorInt(std::ostream &os, std::vector<int> const& OneInc);

struct PairEltPerm {
  permlib::dom_int eElt;
  permlib::Permutation ePerm;
};


struct MyFormTransversal {
  permlib::dom_int eElt;
  std::vector<PairEltPerm> ListPair;
};


permlib::Permutation IdentityPermutation(int const& n);


MyFormTransversal GetListPermutation(PermutationGroupPtr TheGRP,
				     permlib::SchreierTreeTransversal<permlib::Permutation> const& eTrans);

struct IteratorGrp {
  int n;
  std::vector<size_t> ListPos;
  std::vector<MyFormTransversal> ListTrans;
};


template<typename Tint>
IteratorGrp GetInitialIterator(TheGroupFormat<Tint> const& eGRP)
{
  int n=eGRP.n_act();
  std::vector<MyFormTransversal> ListTrans;
  for (auto & eTrans : eGRP.get_group()->U) {
    MyFormTransversal TheTrans=GetListPermutation(eGRP.get_group(), eTrans);
    ListTrans.push_back(TheTrans);
  }
  int nbClass=ListTrans.size();
  std::vector<size_t> ListPos(nbClass, 0);
  return {n, ListPos, ListTrans};
}

template<typename Tint>
struct OrbitMinimumArr {
  int n;
  Tint GRPsize;
  std::vector<MyFormTransversal> ListTrans;
};

template<typename Tint>
OrbitMinimumArr<Tint> GetInitialMinimumArray(TheGroupFormat<Tint> const& eGRP)
{
  //  std::cerr << "GetInitialMinimumArray |GRP|=" << eGRP.size << "\n";
  IteratorGrp eIter=GetInitialIterator(eGRP);
  return {eIter.n, eGRP.size(), eIter.ListTrans};
}

template<typename Tint>
struct ResultMinimum {
  Face eMin;
  Tint OrbitSize;
};

template<typename Tint>
ResultMinimum<Tint> GetMinimumRecord(OrbitMinimumArr<Tint> const& ArrMin, Face const& eFace)
{
  int nbTrans=ArrMin.ListTrans.size();
  int n=ArrMin.n;
  std::vector<Face> ListFace(nbTrans+1,eFace);
  std::vector<size_t> ListPos(nbTrans,0);
  Face FaceMin=eFace;
  Tint nbAtt=0;
  //  std::cerr << "|GRP|=" << ArrMin.GRPsize << "\n";
  auto Increment=[&]() -> int {
    for (int i=0; i<nbTrans; i++)
      if (ListPos[i] < ArrMin.ListTrans[i].ListPair.size() -1) {
	ListPos[i]++;
	for (int j=0; j<i; j++)
	  ListPos[j]=0;
	for (int j=i; j>=0; j--) {
	  int ePos=ListPos[j];
	  for (int iCol=0; iCol<n; iCol++) {
	    int jCol=ArrMin.ListTrans[j].ListPair[ePos].ePerm.at(iCol);
	    ListFace[j][jCol]=ListFace[j+1][iCol];
	  }
	}
	return 0;
      }
    return -1;
  };
  while(true) {
    if (ListFace[0] < FaceMin) {
      FaceMin=ListFace[0];
      nbAtt=0;
    }
    if (ListFace[0] == FaceMin)
      nbAtt++;
    int test=Increment();
    if (test == -1)
      break;
  }
  //  std::cerr << "nbAtt=" << nbAtt << "\n";
  Tint eOrbitSize=ArrMin.GRPsize / nbAtt;
  return {FaceMin, eOrbitSize};
}

permlib::Permutation GetPermutation(IteratorGrp const& eIter);
bool IsFinalIterator(IteratorGrp const& eIter);
int IteratorIncrease(IteratorGrp & eIter);
bool TestBelongingInGroup(IteratorGrp const& eIter, permlib::Permutation const& ePerm);


template<typename Tint>
bool IsSubgroup(TheGroupFormat<Tint> const& g1, TheGroupFormat<Tint> const& g2)
{
  IteratorGrp eIter=GetInitialIterator(g1);
  for (auto & eGen : g2.group->S) {
    bool test=TestBelongingInGroup(eIter, *eGen);
    if (!test)
      return false;
  }
  return true;
}

























#endif
