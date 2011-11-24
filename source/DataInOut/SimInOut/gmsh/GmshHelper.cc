// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "GmshHelper.hh"

namespace CoupledField {

  std::map<Elem::FEType, GmshHelper::ElemNodeMapType> GmshHelper::elemNodeMap_;
  
  //*****************
  //   Constructor
  //*****************
  GmshHelper::GmshHelper() {}


  // **********************
  //   Default Destructor
  // **********************
  GmshHelper::~GmshHelper() {}
  
  void GmshHelper::ElemType2FEType(UInt eType, Elem::FEType& feType, UInt& numNodes)
  {
    switch(eType)
    {
    case 1:
      feType = Elem::ET_LINE2;
      break;
    case 2:
      feType = Elem::ET_TRIA3;
      break;
    case 3:
      feType = Elem::ET_QUAD4;
      break;
    case 4:
      feType = Elem::ET_TET4;
      break;
    case 5:
      feType = Elem::ET_HEXA8;
      break;
    case 6:
      feType = Elem::ET_WEDGE6;
      break;
    case 7:
      feType = Elem::ET_PYRA5;
      break;
    case 8:
      feType = Elem::ET_LINE3;
      break;
    case 9:
      feType = Elem::ET_TRIA6;
      break;
    case 10:
      feType = Elem::ET_QUAD9;
      break;
    case 11:
      feType = Elem::ET_TET10;
      break;
    case 12:
      feType = Elem::ET_HEXA27;
      break;
    case 13:
      feType = Elem::ET_UNDEF;
      numNodes = 18;
      break;
    case 14:
      feType = Elem::ET_UNDEF;
      numNodes = 14;
      break;
    case 15:
      feType = Elem::ET_POINT;
      numNodes = 1;
      break;
    case 16:
      feType = Elem::ET_QUAD8;
      break;
    case 17:
      feType = Elem::ET_HEXA20;
      break;
    case 18:
      feType = Elem::ET_WEDGE15;
      break;
    case 19:
      feType = Elem::ET_PYRA13;
      break;
    case 20:
      feType = Elem::ET_UNDEF;
      numNodes = 9;
      break;
    case 21:
      feType = Elem::ET_UNDEF;
      numNodes = 10;
      break;
    case 22:
      feType = Elem::ET_UNDEF;
      numNodes = 12;
      break;
    case 23:
    case 24:
      feType = Elem::ET_UNDEF;
      numNodes = 15;
      break;
    case 25:
      feType = Elem::ET_UNDEF;
      numNodes = 21;
      break;
    case 26:
      feType = Elem::ET_UNDEF;
      numNodes = 4;
    case 27:
      feType = Elem::ET_UNDEF;
      numNodes = 5;
      break;
    case 28:
      feType = Elem::ET_UNDEF;
      numNodes = 6;
      break;
    case 29:
      feType = Elem::ET_UNDEF;
      numNodes = 20;
      break;
    case 30:
      feType = Elem::ET_UNDEF;
      numNodes = 35;
      break;
    case 31:
      feType = Elem::ET_UNDEF;
      numNodes = 56;
      break;
    default:
      EXCEPTION("Gmsh element type " << eType << " not supported!\n"
                << "Please refer to the Gmsh documentation at\n"
                << "http://www.geuz.org/gmsh for available element types.\n"
                << "Documentation for CFS++ element types may be found in\n"
                << "CFS_SOURCE_DIR/share/doc/developer/AddDoc/element_types.");
      break;
    }

    if(feType != Elem::ET_UNDEF) {
      numNodes = Elem::shapes[feType].numNodes;
    } else {
      WARN("Cannot handle Gmsh element type " << eType 
           << " which has " << numNodes << " nodes.");
    }
  }

  void GmshHelper::FEType2ElemType(Elem::FEType feType, UInt& eType, UInt& numNodes)
  {
    switch(feType)
    {
    case Elem::ET_POINT:
      eType = 15;
      break;
    case Elem::ET_LINE2:
      eType = 1;
      break;
    case Elem::ET_TRIA3:
      eType = 2;
      break;
    case Elem::ET_QUAD4:
      eType = 3;
      break;
    case Elem::ET_TET4:
      eType = 4;
      break;
    case Elem::ET_HEXA8:
      eType = 5;
      break;
    case Elem::ET_WEDGE6:
      eType = 6;
      break;
    case Elem::ET_PYRA5:
      eType = 7;
      break;
    case Elem::ET_LINE3:
      eType = 8;
      break;
    case Elem::ET_TRIA6:
      eType = 9;
      break;
    case Elem::ET_QUAD9:
      eType = 10;
      break;
    case Elem::ET_TET10:
      eType = 11;
      break;
    case Elem::ET_HEXA27:
      eType = 12;
      break;
    case Elem::ET_QUAD8:
      eType = 16;
      break;
    case Elem::ET_HEXA20:
      eType = 17;
      break;
    case Elem::ET_WEDGE15:
      eType = 18;
      break;
    case Elem::ET_PYRA13:
      eType = 19;
      break;
    default:
      EXCEPTION("Wrong FEType for Gmsh.");
      break;
    }

    numNodes = Elem::shapes[feType].numNodes;
  }

  void GmshHelper::InitElemNodeMap()
  {
    // Initialize node mappings from/to Gmsh element types
    // Left is Gmsh node, right is CFS++ node.
    
    typedef ElemNodeMapType::value_type emvt;

    elemNodeMap_[Elem::ET_POINT].insert( emvt(0, 0) );

    elemNodeMap_[Elem::ET_LINE2].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_LINE2].insert( emvt(1, 1) );

    elemNodeMap_[Elem::ET_LINE3].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_LINE3].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_LINE3].insert( emvt(2, 2) );

    elemNodeMap_[Elem::ET_TRIA3].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_TRIA3].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_TRIA3].insert( emvt(2, 2) );

    elemNodeMap_[Elem::ET_TRIA6].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_TRIA6].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_TRIA6].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_TRIA6].insert( emvt(3, 3) );
    elemNodeMap_[Elem::ET_TRIA6].insert( emvt(4, 4) );
    elemNodeMap_[Elem::ET_TRIA6].insert( emvt(5, 5) );

    elemNodeMap_[Elem::ET_QUAD4].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_QUAD4].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_QUAD4].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_QUAD4].insert( emvt(3, 3) );

    elemNodeMap_[Elem::ET_QUAD8].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_QUAD8].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_QUAD8].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_QUAD8].insert( emvt(3, 3) );
    elemNodeMap_[Elem::ET_QUAD8].insert( emvt(4, 4) );
    elemNodeMap_[Elem::ET_QUAD8].insert( emvt(5, 5) );
    elemNodeMap_[Elem::ET_QUAD8].insert( emvt(6, 6) );
    elemNodeMap_[Elem::ET_QUAD8].insert( emvt(7, 7) );

    elemNodeMap_[Elem::ET_QUAD9].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_QUAD9].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_QUAD9].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_QUAD9].insert( emvt(3, 3) );
    elemNodeMap_[Elem::ET_QUAD9].insert( emvt(4, 4) );
    elemNodeMap_[Elem::ET_QUAD9].insert( emvt(5, 5) );
    elemNodeMap_[Elem::ET_QUAD9].insert( emvt(6, 6) );
    elemNodeMap_[Elem::ET_QUAD9].insert( emvt(7, 7) );
    elemNodeMap_[Elem::ET_QUAD9].insert( emvt(8, 8) );

    elemNodeMap_[Elem::ET_TET4].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_TET4].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_TET4].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_TET4].insert( emvt(3, 3) );

    elemNodeMap_[Elem::ET_TET10].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_TET10].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_TET10].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_TET10].insert( emvt(3, 3) );
    elemNodeMap_[Elem::ET_TET10].insert( emvt(4, 4) );
    elemNodeMap_[Elem::ET_TET10].insert( emvt(5, 5) );
    elemNodeMap_[Elem::ET_TET10].insert( emvt(6, 6) );
    elemNodeMap_[Elem::ET_TET10].insert( emvt(7, 7) );
    elemNodeMap_[Elem::ET_TET10].insert( emvt(8, 9) );
    elemNodeMap_[Elem::ET_TET10].insert( emvt(9, 8) );

    elemNodeMap_[Elem::ET_HEXA8].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_HEXA8].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_HEXA8].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_HEXA8].insert( emvt(3, 3) );
    elemNodeMap_[Elem::ET_HEXA8].insert( emvt(4, 4) );
    elemNodeMap_[Elem::ET_HEXA8].insert( emvt(5, 5) );
    elemNodeMap_[Elem::ET_HEXA8].insert( emvt(6, 6) );
    elemNodeMap_[Elem::ET_HEXA8].insert( emvt(7, 7) );

    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(3, 3) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(4, 4) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(5, 5) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(6, 6) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(7, 7) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(8, 8) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(9, 11) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(10, 16) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(11, 9) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(12, 17) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(13, 10) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(14, 18) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(15, 19) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(16, 12) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(17, 15) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(18, 13) );
    elemNodeMap_[Elem::ET_HEXA20].insert( emvt(19, 14) );

    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(3, 3) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(4, 4) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(5, 5) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(6, 6) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(7, 7) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(8, 8) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(9, 11) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(10, 16) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(11, 9) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(12, 17) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(13, 10) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(14, 18) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(15, 19) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(16, 12) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(17, 15) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(18, 13) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(19, 14) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(20, 24) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(21, 20) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(22, 23) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(23, 21) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(24, 22) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(25, 25) );
    elemNodeMap_[Elem::ET_HEXA27].insert( emvt(26, 26) );


    elemNodeMap_[Elem::ET_PYRA5].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_PYRA5].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_PYRA5].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_PYRA5].insert( emvt(3, 3) );
    elemNodeMap_[Elem::ET_PYRA5].insert( emvt(4, 4) );

    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(3, 3) );
    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(4, 4) );
    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(5, 5) );
    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(6, 8) );
    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(7, 9) );
    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(8, 6) );
    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(9, 10) );
    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(10, 7) );
    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(11, 11) );
    elemNodeMap_[Elem::ET_PYRA13].insert( emvt(12, 12) );

//     elemNodeMap_[Elem::PYRA13].insert( emvt(0, 0) );
//     elemNodeMap_[Elem::PYRA13].insert( emvt(1, 1) );
//     elemNodeMap_[Elem::PYRA13].insert( emvt(2, 2) );
//     elemNodeMap_[Elem::PYRA13].insert( emvt(3, 3) );
//     elemNodeMap_[Elem::PYRA13].insert( emvt(4, 4) );
//     elemNodeMap_[Elem::PYRA13].insert( emvt(5, 5) );
//     elemNodeMap_[Elem::PYRA13].insert( emvt(6, 6) );
//     elemNodeMap_[Elem::PYRA13].insert( emvt(7, 9) );
//     elemNodeMap_[Elem::PYRA13].insert( emvt(8, 8) );
//     elemNodeMap_[Elem::PYRA13].insert( emvt(9, 10) );
//     elemNodeMap_[Elem::PYRA13].insert( emvt(10, 7) );
//     elemNodeMap_[Elem::PYRA13].insert( emvt(11, 11) );
//     elemNodeMap_[Elem::PYRA13].insert( emvt(12, 12) );

    elemNodeMap_[Elem::ET_WEDGE6].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_WEDGE6].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_WEDGE6].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_WEDGE6].insert( emvt(3, 3) );
    elemNodeMap_[Elem::ET_WEDGE6].insert( emvt(4, 4) );
    elemNodeMap_[Elem::ET_WEDGE6].insert( emvt(5, 5) );

    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(0, 0) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(1, 1) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(2, 2) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(3, 3) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(4, 4) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(5, 5) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(6, 6) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(7, 8) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(8, 12) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(9, 7) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(10, 13) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(11, 14) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(12, 9) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(13, 11) );
    elemNodeMap_[Elem::ET_WEDGE15].insert( emvt(14, 10) );
  }  

}
