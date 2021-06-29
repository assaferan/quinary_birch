#include "GRAPH_GraphicalFunctions.h"

GraphFunctional::GraphFunctional(int const& inpNbVert, std::function<bool(int const&,int const&)> const& eFct) : nbVert(inpNbVert), f(eFct)
{
  HasVertexColor=false;
}

GraphFunctional::~GraphFunctional()
{
}

GraphFunctional::GraphFunctional(GraphFunctional const& eG)
{
  nbVert=eG.GetNbVert();
  f=eG.GetFCT();
  HasVertexColor=eG.GetHasVertexColor();
  fColor=eG.GetFColor();
}

GraphFunctional GraphFunctional::operator=(GraphFunctional const& eG)
{
  nbVert=eG.GetNbVert();
  f=eG.GetFCT();
  HasVertexColor=eG.GetHasVertexColor();
  fColor=eG.GetFColor();
  return *this;
}

// lighter stuff
int GraphFunctional::GetNbVert() const
{
  return nbVert;
}

std::function<bool(int const&,int const&)> GraphFunctional::GetFCT() const
{
  return f;
}

bool GraphFunctional::GetHasVertexColor() const
{
  return HasVertexColor;
}

std::function<int(int const&)> GraphFunctional::GetFColor() const
{
  return fColor;
}

//
void GraphFunctional::SetFColor(std::function<int(int const&)> const& inpFColor)
{
  HasVertexColor=true;
  fColor=inpFColor;
}

std::vector<int> GraphFunctional::Adjacency(int const& iVert) const
{
  std::vector<int> retList;
  for (int jVert=0; jVert<nbVert; jVert++)
    if (f(iVert, jVert))
      retList.push_back(jVert);
  return retList;
}

bool GraphFunctional::IsAdjacent(int const& iVert, int const& jVert) const
{
  return f(iVert, jVert);
}

int GraphFunctional::GetColor(int const& iVert) const
{
  return fColor(iVert);
}

MyMatrix<int> CreateEmbedding(int const& StartPoint, std::vector<std::vector<int> > const& ListEdge, std::vector<std::vector<int> > const& ListLabel)
{
  int nbEdge=ListEdge.size();
  //
  int nbVert=0;
  for (auto & eEdge : ListEdge)
    for (auto & eVal : eEdge)
      if (eVal > nbVert)
	nbVert=eVal;
  nbVert++;
  //
  int nbLabel=0;
  for (auto & eLabel : ListLabel)
    for (auto & eVal : eLabel)
      if (eVal > nbLabel)
	nbLabel=eVal;
  nbLabel++;
  //
  MyMatrix<int> Embedding(nbVert, nbLabel);
  for (int i=0; i<nbLabel; i++)
    Embedding(StartPoint,i)=0;
  std::vector<int> ListStatus(nbVert,0);
  ListStatus[StartPoint]=1;
  while(true) {
    bool IsFinished=true;
    for (int iEdge=0; iEdge<nbEdge; iEdge++) {
      std::vector<int> eEdge=ListEdge[iEdge];
      int eStat1=ListStatus[eEdge[0]];
      int eStat2=ListStatus[eEdge[1]];
      int RelPos=-1;
      int eVert1=-1, eVert2=-1;
      if (eStat1 == 0 && eStat2 == 1) {
	RelPos=0;
	eVert1=eEdge[1];
	eVert2=eEdge[0];
      }
      if (eStat2 == 0 && eStat1 == 1) {
	RelPos=0;
	eVert1=eEdge[0];
	eVert2=eEdge[1];
      }
      if (RelPos != -1) {
	IsFinished=false;
	for (int i=0; i<nbLabel; i++)
	  Embedding(eVert2,i) = Embedding(eVert1,i);
	for (auto & eVal : ListLabel[iEdge]) {
	  Embedding(eVert2,eVal) = 1 - Embedding(eVert1,eVal);
	}
	ListStatus[eVert2]=1;
      }
    }
    if (IsFinished)
      break;
  }
  return Embedding;
}

int L1_distance(std::vector<int> const& V1, std::vector<int> const& V2)
{
  size_t siz=V1.size();
  assert(siz == V2.size());
  int dist=0;
  for (size_t i=0; i<siz; i++)
    if (V1[i] != V2[i])
      dist++;
  return dist;
}
