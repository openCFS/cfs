#ifndef FILE_CFS_ELEM_2002
#define FILE_CFS_ELEM_2002

#include "Utils/StdVector.hh"
#include "Elements/basefe.hh"

namespace CoupledField
{

  //! class we store description of element
struct Elem
{
public:
 
  //! global element number
  Integer elemNum; 
  
  //! pointer to BaseElem. FE-characteristics of element
  BaseFE * ptElem;
  
  //! connection array
  StdVector<Integer> connect;
  
  //! name of subdomain, to which this element is belogned
  std::string namesd;

  //! flag for refinement
  Boolean refinementFlag; 

  //! number of refinement for the element
  Integer refinementNumber; 

  //! overloading operator =
  Elem & operator=(const Elem& t);

  //! calculation of diameter of element
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
  if (connect.GetSize()==1)
    Error("This function is not valid for this dimension",__FILE__,__LINE__);
  
  Point<2> a=ptArrayOfNodes[connect[1]];
  Point<2> b=ptArrayOfNodes[connect[2]];
  Point<2> c=ptArrayOfNodes[connect[3]];

  return std::max(dist(a,b),dist(b,c));
}

} // end of namespace
#endif
