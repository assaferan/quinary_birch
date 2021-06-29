#include "Permlib_specific.h"

int OnPoints(int const& i, permlib::Permutation const& elt)
{
  return elt.at(i);
}

//
// Specific implementation with permlib
//


/*
std::vector<permlib::dom_int> GetBaseGroup(TheGroupFormat const& eGRP)
{
  std::vector<permlib::dom_int> eList;
  for (auto & eTrans : eGRP.group->U) {
    permlib::dom_int eElt=eTrans.element();
    eList.push_back(eElt);
  }
  return eList;
}
*/

std::set<int> GetSetFrom_DB(Face const& eList)
{
  int nb=eList.count();
  std::set<int> eSet;
  int aRow=eList.find_first();
  for (int i=0; i<nb; i++) {
    eSet.insert(aRow);
    aRow=eList.find_next(aRow);
  }
  return eSet;
}

Face Face_EquivSimplification(Face const& eFace)
{
  int siz=eFace.size();
  int cnt=eFace.count();
  if (2*cnt > siz) {
    Face retFace(siz);
    for (int i=0; i<siz; i++)
      retFace[i]=1 - eFace[i];
    return retFace;
  }
  else {
    return eFace;
  }
}

std::pair<bool,permlib::Permutation> PERMLIB_TestEquivalenceGeneralNaked(int const& n, PermutationGroupPtr const& group, Face const& eList1, Face const& eList2, int const& eMethod)
{
  permlib::Permutation::ptr mappingElement;
  int nb1=eList1.count();
  int nb2=eList2.count();
  if (nb1 != nb2)
    return {false, {}};
  //
  Face NewList1=Face_EquivSimplification(eList1);
  Face NewList2=Face_EquivSimplification(eList2);
  std::set<int> eSet1=GetSetFrom_DB(NewList1);
  std::set<int> eSet2=GetSetFrom_DB(NewList2);
  if (eMethod == 0)
    mappingElement = permlib::setImage_classic(*group, eSet1.begin(), eSet1.end(), eSet2.begin(), eSet2.end());
  if (eMethod == 1)
    mappingElement = permlib::setImage_partition(*group, eSet1.begin(), eSet1.end(), eSet2.begin(), eSet2.end());
  if (mappingElement) {
    return {true, *mappingElement};
  }
  return {false, {}};
}

Face PERMLIB_Canonicalization(int const& n, PermutationGroupPtr const& group, Face const& eList)
{
  DsetList eListI(n), eListO(n);
  int siz=eList.count();
  int aRow=eList.find_first();
  for (int i=0; i<siz; i++) {
    eListI[aRow]=1;
    aRow=eList.find_next(aRow);
  }
  eListO=smallestSetImage(*group, eListI);
  Face TheRet;
  for (int i=0; i<n; i++)
    if (eListO[i] == 1)
      TheRet[i]=1;
#ifdef DEBUG_GROUP
  bool test=PERMLIB_TestEquivalence(n, *group, eList, TheRet);
  if (!test) {
    std::cerr << "We have major debugging to do\n";
    throw TerminalException{1};
  }
#endif
  return TheRet;
}



PermutationGroupPtr PERMLIB_GetStabilizer_general(PermutationGroupPtr const& group, Face const& eList, int const& opt)
{
  Face NewList=Face_EquivSimplification(eList);
  std::set<int> eSet=GetSetFrom_DB(NewList);
  if (opt == 0)
    return permlib::setStabilizer_classic(*group, eSet.begin(), eSet.end());
  else
    return permlib::setStabilizer_partition(*group, eSet.begin(), eSet.end());
}

std::size_t std::hash<permlib::Permutation>::operator()(const permlib::Permutation & eElt) const
{
  size_t len=eElt.size();
  std::vector<permlib::dom_int> V(len);
  for (size_t i=0; i<len; i++)
    V[i] = eElt.at(i);
  return std::hash<std::vector<permlib::dom_int>>()(V);
} 

void WriteVectorInt(std::ostream &os, std::vector<int> const& OneInc)
{
  int i, siz;
  siz=OneInc.size();
  for (i=0; i<siz; i++)
    os << " " << OneInc[i];
  os << "\n";
}


permlib::Permutation IdentityPermutation(int const& n)
{
  std::vector<permlib::dom_int> v(n);
  for (int i=0; i<n; i++)
    v[i]=i;
  return permlib::Permutation(v);
}


MyFormTransversal GetListPermutation(PermutationGroupPtr TheGRP,
				     permlib::SchreierTreeTransversal<permlib::Permutation> const& eTrans)
{
  permlib::dom_int n=TheGRP->n;
  permlib::dom_int eElt=eTrans.element();
  std::unordered_set<permlib::Permutation, std::hash<permlib::Permutation>> ListPermWork;
  for (std::shared_ptr<permlib::Permutation> & p : eTrans.GetMtransversal() ) {
    if (p) {
      permlib::Permutation ePerm=*p;
      ListPermWork.insert(ePerm);
    }
  }
  std::unordered_set<permlib::dom_int> SetOrbit;
  std::vector<PairEltPerm> ListPair;
  std::vector<bool> StatusDone;
  std::function<void(permlib::dom_int,permlib::Permutation)> fInsert=[&](permlib::dom_int const& eVal, permlib::Permutation const& ePerm) -> void {
    std::unordered_set<permlib::dom_int>::iterator iter=SetOrbit.find(eVal);
    if (iter == SetOrbit.end()) {
      SetOrbit.insert(eVal);
      PairEltPerm ePair{eVal, ePerm};
      ListPair.push_back(ePair);
      StatusDone.push_back(false);
    }
  };
  permlib::Permutation ePerm=IdentityPermutation(n);
  fInsert(eElt, ePerm);
  while(true) {
    bool IsFinished=true;
    int len=ListPair.size();
    for (int i=0; i<len; i++)
      if (!StatusDone[i]) {
	StatusDone[i]=true;
	IsFinished=false;
	for (auto & fPerm : ListPermWork) {
	  permlib::dom_int fVal=fPerm.at(ListPair[i].eElt);
	  permlib::Permutation eProd=ListPair[i].ePerm*fPerm;
	  permlib::dom_int gVal=eProd.at(eElt);
	  if (gVal != fVal) {
	    std::cerr << "Inconsistency here on the permutation product\n";
	    throw TerminalException{1};
	  }
	  fInsert(fVal, eProd);
	}
      }
    if (IsFinished)
      break;
  }
  return {eElt, ListPair};
}

permlib::Permutation GetPermutation(IteratorGrp const& eIter)
{
  int n=eIter.n;
  permlib::Permutation ePerm=IdentityPermutation(n);
  int nbClass=eIter.ListPos.size();
  for (int i=0; i<nbClass; i++) {
    size_t ePos=eIter.ListPos[i];
    ePerm = eIter.ListTrans[i].ListPair[ePos].ePerm*ePerm;
  }
  return ePerm;
}

bool IsFinalIterator(IteratorGrp const& eIter)
{
  int nbClass=eIter.ListPos.size();
  for (int i=0; i<nbClass; i++) {
    size_t siz=eIter.ListTrans[i].ListPair.size();
    if (eIter.ListPos[i] < siz-1)
      return false;
  }
  return true;
}

int IteratorIncrease(IteratorGrp & eIter)
{
  int nbClass=eIter.ListPos.size();
  for (int i=0; i<nbClass; i++)
    if (eIter.ListPos[i] < eIter.ListTrans[i].ListPair.size() -1) {
      eIter.ListPos[i]++;
      for (int j=0; j<i; j++)
	eIter.ListPos[j]=0;
      return 0;
    }
  return -1;
}


bool TestBelongingInGroup(IteratorGrp const& eIter, permlib::Permutation const& ePerm)
{
  permlib::Permutation eWork=ePerm;
  int nbClass=eIter.ListTrans.size();
  std::function<bool(permlib::dom_int,int)> fUpdate=[&](permlib::dom_int const& eElt,int const& i) -> bool {
    permlib::dom_int eImg=eWork.at(eElt);
    int len=eIter.ListTrans[i].ListPair.size();
    for (int j=0; j<len; j++)
      if (eIter.ListTrans[i].ListPair[j].eElt == eImg) {
	eWork = eWork * (~eIter.ListTrans[i].ListPair[j].ePerm);
	return true;
      }
    return false;
  };
  for (int i=0; i<nbClass; i++) {
    permlib::dom_int eElt=eIter.ListTrans[i].eElt;
    bool test1=fUpdate(eElt, i);
    if (!test1)
      return false;
  }
  return eWork.isIdentity();
}






















