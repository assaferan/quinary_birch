#ifndef __SCHREIER_SIMS_H_
#define __SCHREIER_SIMS_H_

#include <memory>
#include <unordered_set>

#include "HashMap.h"

// based on "The Schreier-Sims algorithm for matrix groups", Henrik Baarnhielm (https://arxiv.org/pdf/math/0410593.pdf)

template <class GrpElt, class GSetElt>
class OrbitTree;

template <class GrpElt, class GSetElt>
class OrbitNode {
public:
  OrbitNode(const GrpElt & g, const GSetElt & x)
    : coset_rep(g), elt(x) {}

  // This one is the coset representative of all the way from the root
  GrpElt coset_rep;
  GSetElt elt;
  std::vector<size_t> children_idxs;
  // This is not really necessary here - remove after debugging
  std::vector<GrpElt> children_labels;
};

template <class GrpElt, class GSetElt>
class OrbitTree {
public:

  // create an empty tree
  OrbitTree() = default;

  OrbitTree(const GSetElt & x)
  {
    GrpElt id;
    OrbitNode<GrpElt, GSetElt> > root(id,x);
    _vec.push_back(root);
    _elt_to_node_idx[x] = 0;
  }

  inline bool empty() const {return _vec.empty();}
  
  void addChild(const GSetElt & p1, const GSetElt & p2, const GrpElt & l)
  {
    assert(_vec[_elt_to_node_idx[p1]].elt == p1);

    const OrbitNode<GrpElt,GSetElt> & old_node = _vec[_elt_to_node_idx[p1]];
    OrbitNode<GrpElt,GSetElt> new_node(l * old_node.coset_rep, p2);

    old_node.children_idxs.push_back(_vec.size());
    old_node.children_labels.push_back(l);
    _vec.push_back(new_node);

    return;
  }

  inline bool contains(const GSetElt & x) const
  {return (_elt_to_node_idx.find(x) != _elt_to_node_idx.end());}
  
  inline GrpElt orbitElement(const GSetElt & x) const {
    assert(this->contains(x));
    return _vec[_elt_to_node_idx[x]].coset_rep;
  }
  
  inline GrpElt getSchreierGenerator(const GSetElt & p, const GrpElt & s)
  {
    GrpElt t1 = this->orbitElement(p);
    GrpElt t2 = this->orbitElement(s*p);
    return t1*s*t2;
  }

  inline size_t size(void) const {return _vec.size(); }
  
protected:

  std::vector< OrbitNode<GrpElt, GSetElt> > _vec;
  std::unordered_map< GSetElt, size_t > _elt_to_node_idx; 

};

template<typename R, size_t n>
VectorInt<R,n> newBasePoint(const Isometry<R,n> & g);

template <class GrpElt, class GSetElt>
class Group
{
public:
  typedef std::unordered_set< GrpElt > GeneratorSet;
  typedef HashMap< GSetElt > EltSeq;

  Group(const GeneratorSet & S)
  {
    EltSeq empty_seq;
    this->computeBSGS(S);
  }

  inline size_t size(void) const
  {
    size_t sz = 1;
    for (size_t tree_idx = 0; tree_idx < _T.size(); tree_idx++)
      sz *= _T[tree_idx].size();

    return sz;
  }
  
  void getPartialBSGS(const GeneratorSet & S, const EltSeq & B)
  {
    for (size_t i = 0; i < B.size(); i++)
      _B.add(B.at(i));

    for (GSetElt s : S) {
      if (!s.isIdentity()) {
	EltSeq s_B;
	for (GSetElt x : _B) {
	  s_B.add(s*x);
	}
	if (s_B == _B) {
	  GSetElt point = newBasePoint(s);
	  _B.add(point);
	}
	_S.insert(s);
	_S.insert(s.inverse());
      }
    }
    return;
  }

  void computeBSGS(const GeneratorSet & S) {
    EltSeq empty_seq;
    this->getPartialBSGS(S, empty_seq);
    
    for (size_t i = base.size(); i > 0; i--) {
      this->doSchreierSims(i); 
    }
    
    return;
  }

  std::pair<GrpElt, size_t> isMember(const GrpElt & g, size_t start_idx)
  {
    size_t n = _B.size();
    GrpElt r = g;
    for (size_t i = start_idx; i < n; i++) {
      GSetElt x = r * _B[i];
      if (!(_T[i].contains(x))) {
	return make_pair(r, i);
      }
      assert(i < _T.size());
      GrpElt elt = _T[i].orbitElement(x);
      r = r*elt.inverse();
    }
    return make_pair(r,n);
  }
  
  void doSchreierSims(size_t i) {
    GeneratorSet gens(_S.begin(), _S.end()); // Is this S_i ?
    size_t k = _B.size();
    assert( i <= k );
    if (_T.size() < k)
      _T.resize(k);
    OrbitTree<GrpElt, GSetElt> T = computeSchreierTree(gens, _B[i-1]);
    for (GSetElt p : T.orbit()) {
      for (GrpElt s : gens) {
	GrpElt gen = T.getSchreierGenerator(p,s);
	if (!gen.isIdentity()) {
	  std::pair<GrpElt, size_t> is_member = this->isMember(gen, i);
	  GrpElt residue = is_member.first;
	  size_t dropout = is_member.second;
	  if (!residue.isIdentity()) {
	    _S.insert(gen);
	    _S.insert(gen.inverse());
	    if (dropout == k) {
	      GSetElt point = newBasePoint(gen);
	      _B.insert(point);
	    }
	    for (size_t j = dropout; j >= i; j--) {
	      doSchreierSims(j+1);
	    }
	  }
	}
      }
    }
    _T[i-1] = T;
    return;
  }

  OrbitTree<GrpElt, GSetElt> computeSchreierTree(const GeneratorSet & S, const GSetElt & x)
  {
    OrbitTree<GrpElt, GSetElt> tree(x);
    EltSeq points;
    points.add(x);

    do {
      EltSeq children;
      for (size_t idx = 0; idx < points.size(); idx++) {
	GSetElt p = points.at(idx);
	for (GrpElt s : S) {
	  GSetElt p_prime = s * p;
	  if (!tree.contains(p_prime)) {
	    tree.addChild(p, p_prime, s);
	    children.add(p_prime);
	  }
	}
      }
      points = children;
    } while (points.size() > 0);

    return tree;
  }
  
  
protected:
  GeneratorSet _S;
  EltSeq _B;
  std::vector< OrbitTree<GrpElt, GSetElt> > _T;
};


template<typename R, size_t n>
VectorRat<R,n> newBasePoint(const Isometry<R,n> & g)
{
  std::shared_ptr< const RationalField<R> > QQ = std::make_shared< const RationalField<R> >();
  VectorRat<R,n> v(QQ);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      if ((i != j) && (g(i,j) != 0)) {
	v[i].makeOne();
	assert(g*v != v);
	return v;
      }

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      if ((i != j) && (g(i,i) != g(j,j))) {
	v[i].makeOne();
	v[j].makeOne();
	assert(g*v != v);
	return v;
      }
  
  v[0].makeOne();
  assert(g*v != v);
  return v;
}

/*
template<typename R, typename S, typename T, size_t n>
VectorFp<R,S,n> Group< Isometry<T,n>, VectorFp<R,S,n> >::newBasePoint(const Isometry<T,n> & g, std::shared_ptr< Fp<R,S> > GF)
{
  VectorFp<R,S,n> v(GF);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      if ((i != j) && (g(i,j) != 0)) {
	v[i].makeOne();
	assert(g*v != v);
	return v;
      }

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      if ((i != j) && (g(i,i) != g(j,j))) {
	v[i].makeOne();
	v[j].makeOne();
	assert(g*v != v);
	return v;
      }
  
  v[0].makeOne();
  assert(g*v != v);
  return v;
}
*/

#endif //__SCHREIER_SIMS_H_
