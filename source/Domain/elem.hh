// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_ELEM_2002
#define FILE_CFS_ELEM_2002

#include <bitset>

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
    Elem() : 
      elemNum(0),
      regionId( NO_REGION_ID )
      {;}

    //! Dummy destructor
    virtual ~Elem() {;}

    // ======================================================
    // GEOMETRICAL INFORMATION
    // ======================================================

    //@{ \name Geometrical Information
    //! global element number
    UInt elemNum; 

    //! identifier for region
    RegionIdType regionId;

    //! array with node numbers
    StdVector<UInt> connect;

    //! array with edge numbers
    StdVector<Integer> edges;

    //! array with face numbers
    StdVector<Integer> faces;

    //! bitset describing the orientation of the faces (3 for each)
    StdVector<std::bitset<3> > faceFlags;
  
#ifdef ADAPTGRID
    //! flag for refinement
    bool refinementFlag; 
  
    //! number of refinement for the element
    UInt refinementNumber; 

#endif
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
    Double diameter(const Point * const ptArrayOfNodes);
    //@}
    
    std::string ToString() const
    {
      std::ostringstream os;
      os << "elemNum=" << elemNum << " region=" << regionId;
      return os.str();
    }
  };



  inline Elem & Elem::operator=(const Elem& t) 
  {
    if (this!=&t) {
      ptElem=t.ptElem;
      connect=t.connect;
      regionId=t.regionId;
#ifdef ADAPTGRID
      refinementFlag=t.refinementFlag;
      refinementNumber=t.refinementNumber;
#endif
    }
    return *this;
  }

  inline Double Elem::diameter(const Point * const ptArrayOfNodes)
  {
    if (connect.GetSize()==1)
      Error("This function is not valid for this dimension",__FILE__,__LINE__);
  
    Point a=ptArrayOfNodes[connect[1]];
    Point b=ptArrayOfNodes[connect[2]];
    Point c=ptArrayOfNodes[connect[3]];

    return std::max(Point::dist(a,b), Point::dist(b,c));
  }

} // end of namespace
#endif
