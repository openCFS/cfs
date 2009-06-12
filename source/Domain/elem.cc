// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "elem.hh"
#include "Elements/basefe.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "General/Enum.hh"

namespace CoupledField {

  // Definition of finite element type mappings
  static EnumTuple feTypeTuples[] = {
    EnumTuple(Elem::UNDEF,   "UNDEF"),
    EnumTuple(Elem::POINT,   "POINT"),
    EnumTuple(Elem::LINE2,   "LINE2"),
    EnumTuple(Elem::LINE3,   "LINE3"),
    EnumTuple(Elem::TRIA3,   "TRIA3"),
    EnumTuple(Elem::TRIA6,   "TRIA6"),
    EnumTuple(Elem::QUAD4,   "QUAD4"),
    EnumTuple(Elem::QUAD8,   "QUAD8"),
    EnumTuple(Elem::QUAD9,   "QUAD9"),
    EnumTuple(Elem::TET4,    "TET4"),
    EnumTuple(Elem::TET10,   "TET10"),
    EnumTuple(Elem::HEXA8,   "HEXA8"),
    EnumTuple(Elem::HEXA20,  "HEXA20"),
    EnumTuple(Elem::HEXA27,  "HEXA27"),
    EnumTuple(Elem::PYRA5,   "PYRA5"),
    EnumTuple(Elem::PYRA13,  "PYRA13"),
    EnumTuple(Elem::WEDGE6,  "WEDGE6"),
    EnumTuple(Elem::WEDGE15, "WEDGE15")      
  };

  Enum<Elem::FEType> Elem::feType = \
  Enum<Elem::FEType>("Finite Element Geometry Types",
                     sizeof(feTypeTuples) / sizeof(EnumTuple),
                     feTypeTuples); 

  
  std::map<Elem::FEType, UInt> Elem::numElemNodes_;
  std::map<Elem::FEType, UInt> Elem::elemDims_;
  std::map<Elem::FEType, UInt> Elem::elemQuadratic_;

  Elem::Elem() : 
    elemNum(0),
    regionId( NO_REGION_ID ),
    ptElem(NULL)
  {
    Initialize();
  }

  void Elem::Initialize()
  {
    if(!numElemNodes_.empty())
      return;

    numElemNodes_[UNDEF] =  0;
    numElemNodes_[POINT] =  1;
    numElemNodes_[LINE2] =  2;
    numElemNodes_[LINE3] =  3;
    numElemNodes_[TRIA3] =  3;
    numElemNodes_[TRIA6] =  6;
    numElemNodes_[QUAD4] =  4;
    numElemNodes_[QUAD8] =  8;
    numElemNodes_[QUAD9] =  9;
    numElemNodes_[TET4]  =  4;
    numElemNodes_[TET10] =  10;
    numElemNodes_[HEXA8] =  8;
    numElemNodes_[HEXA20]=  20;
    numElemNodes_[HEXA27]=  27;
    numElemNodes_[PYRA5] =  5;
    numElemNodes_[PYRA13]=  13;
    numElemNodes_[WEDGE6]=  6;
    numElemNodes_[WEDGE15]= 15;      

    elemDims_[UNDEF] =   0;
    elemDims_[POINT] =   0;
    elemDims_[LINE2] =   1;
    elemDims_[LINE3] =   1;
    elemDims_[TRIA3] =   2;
    elemDims_[TRIA6] =   2;
    elemDims_[QUAD4] =   2;
    elemDims_[QUAD8] =   2;
    elemDims_[QUAD9] =   2;
    elemDims_[TET4] =    3;
    elemDims_[TET10] =   3;
    elemDims_[HEXA8] =   3;
    elemDims_[HEXA20] =  3;
    elemDims_[HEXA27] =  3;
    elemDims_[PYRA5] =   3;
    elemDims_[PYRA13] =  3;
    elemDims_[WEDGE6] =  3;
    elemDims_[WEDGE15] = 3;      

    elemQuadratic_[UNDEF] =   false;
    elemQuadratic_[POINT] =   false;
    elemQuadratic_[LINE2] =   false;
    elemQuadratic_[LINE3] =   true;
    elemQuadratic_[TRIA3] =   false;
    elemQuadratic_[TRIA6] =   true;
    elemQuadratic_[QUAD4] =   false;
    elemQuadratic_[QUAD8] =   true;
    elemQuadratic_[QUAD9] =   true;
    elemQuadratic_[TET4] =    false;
    elemQuadratic_[TET10] =   true;
    elemQuadratic_[HEXA8] =   false;
    elemQuadratic_[HEXA20] =  true;
    elemQuadratic_[HEXA27] =  true;
    elemQuadratic_[PYRA5] =   false;
    elemQuadratic_[PYRA13] =  true;
    elemQuadratic_[WEDGE6] =  false;
    elemQuadratic_[WEDGE15] = true;   
  }

  UInt Elem::GetNumElemNodes(Elem::FEType type) {
    if(numElemNodes_.empty())
      Initialize();

    if(numElemNodes_.find(type) != numElemNodes_.end())
      return numElemNodes_[type];
    else
      return 0;
  } 

  
  UInt Elem::GetElemDim(Elem::FEType type) {
    if(numElemNodes_.empty())
      Initialize();

    if(elemDims_.find(type) != elemDims_.end())
      return elemDims_[type];
    else
      return 0;
  }

  bool Elem::GetElemQuadratic(FEType type) 
  {
    if(numElemNodes_.empty())
      Initialize();

    if(elemQuadratic_.find(type) != elemQuadratic_.end())
      return elemQuadratic_[type];
    else
      return 0;
  }

  void Elem::CorrectConnectivity() 
  {
    UInt dummy;
    FEType type = ptElem->feType();

    // This funtion rearanges the connectivity of the element so
    // that the orientation is in such a way that the Jacobion determinant
    // will always be positive. 
    switch(type) 
    {
    case TRIA6:
      dummy = connect[4];
      connect[4] = connect[5];
      connect[5] = dummy;
    case TRIA3:
      dummy = connect[0];
      connect[0] = connect[1];
      connect[1] = dummy;
      break;
    case QUAD8:
      dummy = connect[4];
      connect[4] = connect[7];
      connect[7] = dummy;
      dummy = connect[5];
      connect[5] = connect[6];
      connect[6] = dummy;
    case QUAD4:
      dummy = connect[1];
      connect[1] = connect[3];
      connect[3] = dummy;
      break;
    default:
      EXCEPTION("Connectivity for " << feType.ToString(type) << " element "
                << elemNum << " in region "
                << domain->GetGrid()->RegionIdToName(regionId)
                << " is not properly oriented!" );
      break;
    }
  }
  
}
