#ifndef FILE_CFS_ELEM_2002
#define FILE_CFS_ELEM_2002

#include "baseelem.hh"

namespace CoupledField
{

struct Elem
{
  BaseElem * ptElem;
  Vector<Integer> connect;
  std::string namesd;
  //! flag for refinement; TRUE - then this element should be refined
  Boolean refinementFlag;

  Elem & operator=(const Elem& t);
  Double diameter(const Point<2> * const ptArrayOfNodes);
};

inline Elem & Elem::operator=(const Elem& t) 
{
  if (this!=&t) {
    ptElem=t.ptElem;
    connect=t.connect;
    namesd=t.namesd;
  }
  return *this;
}

inline Double Elem::diameter(const Point<2> * const ptArrayOfNodes)
{
  if (connect.size()==1)
    Error("This function is not valid for this dimension",__FILE__,__LINE__);
  
  Point<2> a=ptArrayOfNodes[connect[1]];
  Point<2> b=ptArrayOfNodes[connect[2]];
  Point<2> c=ptArrayOfNodes[connect[3]];
    
  return std::max(dist(a,b),dist(b,c));
}

} // end of namespace
#endif
