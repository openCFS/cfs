#ifndef FILE_CFS_ELEM_2002
#define FILE_CFS_ELEM_2002

#include "Utils/StdVector.hh"
#include "Elements/basefe.hh"

namespace CoupledField
{

  //! Class for description of a volume finite element

  //! This class describes a volume finite element, where volume means the 
  //! highest dimensional element entities in the current mesh.
  //! It has to be very lightweight, since this object is created many times.
  //! It relates the geometric information of an element (node numbers)
  //! with the mathematical / computational one (reference finite element).
  //! The finite element is described by:
  //! - corner node numbers
  //! - pointer to reference finite element
  //! - element number
  //! - element subdomain identifier
  //! - refinement flag / number
  

  struct Elem
  {
  public:
 
    //! Dummy constructor
    Elem() {;}

    //! Dummy destructor
    virtual ~Elem() {;}

    // ======================================================
    // GEOMETRICAL INFORMATION
    // ======================================================

    //@{ \name Geometrical Information
    //! global element number
    Integer elemNum; 

    //! identifier for region
    Integer regionId;

    //! array with node numbers
    StdVector<Integer> connect;
  
    //! flag for refinement
    Boolean refinementFlag; 
  
    //! number of refinement for the element
    Integer refinementNumber; 
  
    //@}

    // ======================================================
    // COMPUTATIONAl INFORMATION
    // ======================================================  
  
    //@{ \name Computational Information

    //! pointer to reference element representation
    BaseFE * ptElem;
    //@}

    // ======================================================
    // HELPER METHODS
    // ======================================================
    //@{ \name Helper Methods
  
    //! overloading operator =
    Elem & operator=(const Elem& t);

    //! calculation of diameter of element
    Double diameter(const Point<2> * const ptArrayOfNodes);
    //@}
  };



  inline Elem & Elem::operator=(const Elem& t) 
  {
    if (this!=&t) {
      ptElem=t.ptElem;
      connect=t.connect;
      regionId=t.regionId;
      refinementFlag=t.refinementFlag;
      refinementNumber=t.refinementNumber;
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
