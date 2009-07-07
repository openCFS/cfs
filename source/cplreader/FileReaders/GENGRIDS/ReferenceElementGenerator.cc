// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <cmath>
#include <iterator>

#include "Domain/elem.hh"
#include "cplreader/Settings.hh"
#include "ReferenceElementGenerator.hh"

namespace CoupledField
{
  ReferenceElementGenerator::ReferenceElementGenerator()
  {
  }
  
  void ReferenceElementGenerator::GenGrid(std::vector<double>& coords,
                                          std::vector<UInt>& connect,
                                          std::vector<UInt>& elemTypes,
                                          UInt& maxNumElemNodes,
                                          std::map<std::string, std::vector<UInt> >& regionElems,
                                          std::map<std::string, std::vector<UInt> >& nodeGroups,
                                          Elem::FEType elemType)
  {
    coords.clear();
    connect.clear();

    maxNumElemNodes = Elem::GetNumElemNodes(elemType);
    UInt idx;
    
    elemTypes.clear();

    regionElems.clear();
    std::vector<UInt> elems;
    elems.push_back(1);
    regionElems[Elem::feType.ToString(elemType)] = elems;
    elemTypes.resize(1);
    elemTypes[0] = elemType;
    coords.resize(maxNumElemNodes * 3);
    connect.resize(maxNumElemNodes);
    
    switch(elemType) 
    {
    case Elem::LINE3:
      idx = 2;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
    case Elem::LINE2:
      idx = 0;
      coords[idx*3+0] = -1.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 1;
      coords[idx*3+0] = 1.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      break;

    case Elem::TRIA6:
      idx = 3;
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 4;
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 5;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
    case Elem::TRIA3:
      idx = 0;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 1;
      coords[idx*3+0] = 1.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 2;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = 1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      break;

    case Elem::QUAD9:
      idx = 8;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
    case Elem::QUAD8:
      idx = 4;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = -1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 5;
      coords[idx*3+0] = 1.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 6;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = 1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 7;
      coords[idx*3+0] = -1.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
    case Elem::QUAD4:
      idx = 0;
      coords[idx*3+0] = -1.0;
      coords[idx*3+1] = -1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 1;
      coords[idx*3+0] = 1.0;
      coords[idx*3+1] = -1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 2;
      coords[idx*3+0] = 1.0;
      coords[idx*3+1] = 1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 3;
      coords[idx*3+0] = -1.0;
      coords[idx*3+1] = 1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      break;

    case Elem::TET10:
      idx = 4;
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

      idx = 5;
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

      idx = 6;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

      idx = 7;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;

      idx = 8;
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;

      idx = 9;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;

    case Elem::TET4:
      idx = 0;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

      idx = 1;
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

      idx = 2;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

      idx = 3;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
      break;

    case Elem::HEXA27:
      idx = 20;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 21;          
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 22;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 23;          
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 24;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                          
      idx = 25;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                          
      idx = 26;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

    case Elem::HEXA20:
      idx = 8;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                          
      idx = 9;          
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                          
      idx = 10;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                          
      idx = 11;          
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                          
      idx = 12;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                          
      idx = 13;          
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                          
      idx = 14;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                          
      idx = 15;          
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                          
      idx = 16;          
      coords[idx*3+0] = -1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 17;          
      coords[idx*3+0] = 1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 18;          
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 19;          
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

    case Elem::HEXA8:
      idx = 0;
      coords[idx*3+0] = -1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 1;            
      coords[idx*3+0] = 1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 2;           
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 3;           
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 4;           
      coords[idx*3+0] = -1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;    
                            
      idx = 5;           
      coords[idx*3+0] = 1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;    
                            
      idx = 6;           
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;    
                            
      idx = 7;           
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
      break;

    case Elem::PYRA13:
      idx = 5;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 6;           
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 7;           
      coords[idx*3+0] = 0;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 8;           
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 9;           
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;
                            
      idx = 10;            
      coords[idx*3+0] = -0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;
                             
      idx = 11;            
      coords[idx*3+0] = -0.5;
      coords[idx*3+1] = -0.5;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;
                             
      idx = 12;            
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = -0.5;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;

    case Elem::PYRA5:
      idx = 0;
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 1;            
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 2;           
      coords[idx*3+0] = -1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 3;           
      coords[idx*3+0] = 1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 4;           
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
      break;

    case Elem::WEDGE15:
      idx = 6;           
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                            
      idx = 7;             
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                           
      idx = 8;             
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                           
      idx = 9;             
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                           
      idx = 10;            
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                           
      idx = 11;            
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                           
      idx = 12;            
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                           
      idx = 13;            
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                           
      idx = 14;            
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

    case Elem::WEDGE6:
      idx = 0;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 1;            
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 2;            
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 3;            
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;    
                            
      idx = 4;            
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;    
                            
      idx = 5;            
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;    
      break;

    default:
      break;
    }
    
    
  }
  
} // end of namespace
