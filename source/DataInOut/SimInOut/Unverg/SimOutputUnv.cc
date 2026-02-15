// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <exception>
#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <utility>
#include "SimOutputUnv.hh"
#include "Utils/Timer.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "def_cfs_stats.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  SimOutputUnv::SimOutputUnv(  const std::string& filename,
                               PtrParamNode outputNode,
                               PtrParamNode infoNode, 
                               bool isRestart ) 
    : SimOutput ( filename, outputNode, infoNode, isRestart ),
      output(NULL),
      capaOut_(false) 
  {

    // The restart case is currently not implemented, i.e. resuls from a 
    // partial simulation get overwritten.
    if( isRestart_ ) {
      WARN( "The UNV-Writer is currently not adapted to write restarted "
            "results correctly, thus the results of the previous run get"
            " overwritten." );
    }
    
    
    std::string flavor;
    outputNode->GetValue("flavor", flavor, ParamNode::PASS );
    if(flavor == "CAPA") {
      capaOut_ = true;
    }

    if(capaOut_) {
      formatName_ = "unverg";
    } else {
      formatName_ = "unv";
    }
    
    std::string dirString = "results_" + formatName_; 
    outputNode->GetValue("directory", dirString, ParamNode::PASS );
    dirName_ = dirString; 

    fileName_ =  filename;
    
    stepNumOffset_ = 0;
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );
   
  }


  // **************
  //   Destructor
  // **************
  SimOutputUnv::~SimOutputUnv() {
    delete output;
  }


  // *************
  //   WriteGrid
  // *************
  void SimOutputUnv::WriteGrid() {


    if ( !output ) {
      EXCEPTION( "File for output results is not initialized" );
    }

    if(capaOut_) 
    {
      Dataset666();
      Dataset781();
      Dataset780();
    } else 
    {    
      Dataset151();
      Dataset164();
      Dataset2411();
      Dataset2412();
      Dataset2467();
    }    
  }

  void SimOutputUnv::Dataset151() 
  {
    if (!ptGrid_) {
      EXCEPTION("ptGrid_ is not initialized" );
    }
    char buf[128];
    std::string flavor = capaOut_ == true ? "CAPA" : "SDRC I-DEAS";

    (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 151 << std::endl ;

#if _MSC_VER >= 1400
#define snprintf _snprintf
#endif

    // Model name
    snprintf(buf, 81, "%s", fileName_.c_str());
    (*output) << buf << std::endl;

    // Model description
    snprintf(buf, 81, "UNV Flavor: %s", flavor.c_str());
    (*output) << buf << std::endl;

    // DB Application
    (*output) << "openCFS (TU Wien, TUGRAZ and FAU)" << std::endl;

    std::stringstream sstr;
    sstr << Timer::TimeStamp();
    
    (*output) << sstr.str();
    (*output) << std::setw(10) << 1 
              << std::setw(10) << 0 
              << std::setw(10) << 0
              << std::endl;

    (*output) << sstr.str() << std::endl;

    // UNV Creating Application
    (*output) << "openCFS" << std::endl;

    UInt year, month;
    (*output) << sstr.str();
    sstr.clear(); sstr.str("");
    sstr << CFS_VERSION_YEAR;
    sstr >> year;    
    sstr.clear(); sstr.str("");
    sstr << CFS_VERSION_MONTH;
    sstr >> month;    
    (*output) << std::setw(10) << (year*100 + month);
    sstr.clear(); sstr.str("");
    sstr << CFS_GIT_COMMIT;
    UInt rev;
    std::string modif;
    sstr >> rev >> modif;
    (*output) << std::setw(10) << rev;
    (*output) << std::setw(10) << 0;
    (*output) << std::setw(10) << (modif == "modified" ? 1 : 0);
    (*output) << std::setw(10) << 0 << std::endl;
    
    (*output) << std::setw(6) << -1 << std::endl;
    
  }
  
  void SimOutputUnv::Dataset164()
  {
    // We use SI units in openCFS
    (*output) << std::setw(6) << -1 << std::endl << std::setw(6)
              << 164 << std::endl ;
    (*output) << "         1Meter (newton)               2" << std::endl;
    (*output) << "  1.00000000000000000D+00  1.00000000000000000D+00  1.00000000000000000D+00" << std::endl;
    (*output) << "  2.73149999999999980D+02" << std::endl;
    (*output) << std::setw(6) << -1 << std::endl;
  }

  void  SimOutputUnv::Dataset666() {

    if (!ptGrid_) {
      EXCEPTION("ptGrid_ is not initialized" );
    }

    (*output) << std::setw(6) << -1 << std::endl << std::setw(6)
              << -666 << std::endl ;

    UInt dim=ptGrid_->GetDim();
    UInt maxnumnodes=ptGrid_->GetNumNodes();
    UInt maxnumelem=ptGrid_->GetNumVolElems();

    (*output) << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10)
              << dim << std::endl << std::setw(10) << maxnumnodes
              << std::setw(10) << maxnumelem << std::endl;

    (*output) << std::setw(6) << -1 << std::endl;
  }

  void  SimOutputUnv::Dataset781() {
    //
    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized" );

    (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 781
              << std::endl;

    UInt dim=ptGrid_->GetDim();
    UInt maxnumnodes=ptGrid_->GetNumNodes();

    (*output).setf(std::ios::scientific);
    (*output).precision(16);

    for ( UInt i = 0; i < maxnumnodes; i++ ) {
      (*output) << std::setw(10) << i+1 << std::setw(10) << 0
                << std::setw(10) << 0 << std::setw(10) << 11 << std::endl;

      (*output).setf(std::ios::uppercase);

      Vector<Double> Point;
      ptGrid_->GetNodeCoordinate3D(Point,i+1);
      if( dim == 3 ) {
        (*output) << "   " << Point[0];
        (*output) << "   " << Point[1];
        (*output) << "   " << Point[2];
      } else {
        (*output) << "   " << 0.0;
        (*output) << "   " << Point[0];
        (*output) << "   " << Point[1];
      }
      (*output) << std::endl;
    }
 
    (*output) << std::setw(6) << -1 << std::endl;
  }

  void  SimOutputUnv::Dataset2411() {
    //
    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized" );

    (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 2411
              << std::endl;

    UInt maxnumnodes=ptGrid_->GetNumNodes();

    (*output).setf(std::ios::scientific);
    (*output).precision(16);

    for ( UInt i = 0; i < maxnumnodes; i++ ) {
      (*output) << std::setw(10) << i+1 << std::setw(10) << 1
                << std::setw(10) << 1 << std::setw(10) << 11 << std::endl;

      (*output).setf(std::ios::uppercase);

      Vector<Double> Point;
      Double x,y,z;
      ptGrid_->GetNodeCoordinate3D(Point,i+1);
      x = Point[0];
      y = Point[1];
      z = Point[2];

      // hack to print the numbers with "D+XX" exponents
      char tmp[128];
      sprintf(tmp, "%25.16E%25.16E%25.16E", x, y, z);
      for(unsigned int i = 0; i < strlen(tmp); i++) if(tmp[i] == 'E') tmp[i] = 'D';
      (*output) << tmp << std::endl;
    }
 
    (*output) << std::setw(6) << -1 << std::endl;
  }

  void  SimOutputUnv::Dataset780()
  {
    //
    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized");

    (*output) << std::setw(6) << -1 << std::endl << std::setw(6)
              << 780 << std::endl;

    UInt dim = ptGrid_->GetDim();

    StdVector<UInt> connect;
    StdVector<Elem*> elemssd;
    UInt elmsgrp=1;
    std::string errMsg;

    StdVector<RegionIdType> subdoms;
    ptGrid_->GetVolRegionIds(subdoms);

    for (unsigned int i=0; i<subdoms.GetSize(); i++) {

      ptGrid_->GetVolElems(elemssd,subdoms[i]);

      for (unsigned int j = 0; j < elemssd.GetSize(); j++ ) {
        connect=elemssd[j]->connect;
        (*output) << std::setw(10) << elemssd[j]->elemNum << std::setw(10);

        if ( dim == 2 ) {
          switch(connect.GetSize()) {
          case 3: (*output) << 91; break;
          case 4: (*output) << 94; break;
          case 6: (*output) << 92; break;
          case 8: (*output) << 95; break;
          default:
            EXCEPTION( "Please, put element type according to "
                       << "unverg-format for this number of nodes "
                       << "per element" );
          }

          (*output) << std::setw(10) << 2 << std::setw(10) << 2
                    << std::setw(10) << 1 << std::setw(10) << 1
                    << std::setw(10) << elmsgrp << std::setw(10)
                    << connect.GetSize() << std::endl;
        }
        else {
          switch(connect.GetSize()) {
          case  4: (*output) << 111; break;  // tetraeder 1.ord
          case  6: (*output) << 112; break;  // prism     1.ord
          case  8: (*output) << 115; break;  // hexaeder  1.ord
          case 15: (*output) << 113; break;  // prism     2.ord
          case 20: (*output) << 116; break;  // hexaeder  2.ord
          default:
            EXCEPTION( "Please, put element type according to "
                       << "unverg-format for " << connect.GetSize()
                       << " nodes per Element in 3D!" );
          }

          (*output) << std::setw(10) << 11 << std::setw(10) << 1
                    << std::setw(10) << 1  << std::setw(10) << 1
                    << std::setw(10) << elmsgrp << std::setw(10) 
                    << connect.GetSize() << std::endl;
        }

        if (dim == 2 && (connect.GetSize() == 6 || connect.GetSize() == 8)) {
          //quadratic elements
          UInt offset = UInt(connect.GetSize()/2);
          for (UInt ii=0; ii < offset; ii++) {
            (*output).width(10);
            (*output) << connect[ii];
            (*output).width(10);
            (*output) << connect[ii+offset]; 
          }
        }
        else {
          for (UInt ii=0; ii < connect.GetSize(); ii++) { 
            (*output).width(10);
            (*output) << connect[ii];
          }
        }

        (*output) << std::endl;
      } // over elements of group
      elmsgrp++;
    } // over groups
    (*output) << std::setw(6) << -1 << std::endl;
  }

  void  SimOutputUnv::Dataset2412()
  {
    //
    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized");

    (*output) << std::setw(6) << -1 << std::endl << std::setw(6)
              << 2412 << std::endl;

    StdVector<UInt> connect;
    StdVector<Elem*> elemssd;
    UInt elmsgrp=1;
    std::string errMsg;

    StdVector<RegionIdType> subdoms;
    ptGrid_->GetVolRegionIds(subdoms);
    StdVector<RegionIdType> surfsubdoms;
    ptGrid_->GetSurfRegionIds(surfsubdoms);

    UInt i, j;
    for (i=0; i<surfsubdoms.GetSize(); i++) {
      subdoms.Push_back(surfsubdoms[i]);
    }

    for (i=0; i<subdoms.GetSize(); i++) {

      ptGrid_->GetElems(elemssd,subdoms[i]);

      for ( j = 0; j < elemssd.GetSize(); j++ ) {  
        connect=elemssd[j]->connect;
        UInt unvId = 0;
        
        Elem::FEType elemType = elemssd[j]->type;
        ElemType2UnvElemId(elemType, unvId);

        (*output) << std::setw(10) << elemssd[j]->elemNum;
        (*output) << std::setw(10) << unvId;
        (*output) << std::setw(10) << elmsgrp
                  << std::setw(10) << elmsgrp 
                  << std::setw(10) << elmsgrp // ((UInt) elemType)
                  << std::setw(10);

        switch(elemType) 
        {
        case Elem::ET_LINE2:
          (*output) << std::setw(10) << 2 << std::endl;
          (*output) << std::setw(10) << 0;
          (*output) << std::setw(10) << 1;
          (*output) << std::setw(10) << 1 << std::endl;;
          (*output) << std::setw(10) << connect[0];
          (*output) << std::setw(10) << connect[1];
          break;

        case Elem::ET_LINE3:
          (*output) << std::setw(10) << 2 << std::endl;
          (*output) << std::setw(10) << 0;
          (*output) << std::setw(10) << 1;
          (*output) << std::setw(10) << 1 << std::endl;;
          (*output) << std::setw(10) << connect[0];
          (*output) << std::setw(10) << connect[1];
          (*output) << std::setw(10) << connect[2];
          break;

        case Elem::ET_TRIA6:
          (*output) << std::setw(10) << 6 << std::endl;
          (*output) << std::setw(10) << connect[0];
          (*output) << std::setw(10) << connect[3];
          (*output) << std::setw(10) << connect[1];
          (*output) << std::setw(10) << connect[4];
          (*output) << std::setw(10) << connect[2];
          (*output) << std::setw(10) << connect[5];
          break;

        case Elem::ET_QUAD8:
        case Elem::ET_QUAD9:
          (*output) << std::setw(10) << 8 << std::endl;
          (*output) << std::setw(10) << connect[0];
          (*output) << std::setw(10) << connect[4];
          (*output) << std::setw(10) << connect[1];
          (*output) << std::setw(10) << connect[5];
          (*output) << std::setw(10) << connect[2];
          (*output) << std::setw(10) << connect[6];
          (*output) << std::setw(10) << connect[3];
          (*output) << std::setw(10) << connect[7];
          break;
          
        case Elem::ET_TET10:
          (*output) << std::setw(10) << 10 << std::endl;
          (*output) << std::setw(10) << connect[0];
          (*output) << std::setw(10) << connect[4];

          (*output) << std::setw(10) << connect[1];
          (*output) << std::setw(10) << connect[5];

          (*output) << std::setw(10) << connect[2];
          (*output) << std::setw(10) << connect[6];

          (*output) << std::setw(10) << connect[7];
          (*output) << std::setw(10) << connect[8] << std::endl;
          (*output) << std::setw(10) << connect[9];

          (*output) << std::setw(10) << connect[3];
          break;

        case Elem::ET_PYRA5:
          (*output) << std::setw(10) << 8 << std::endl;
          (*output) << std::setw(10) << connect[0];
          (*output) << std::setw(10) << connect[1];
          (*output) << std::setw(10) << connect[2];
          (*output) << std::setw(10) << connect[3];
          (*output) << std::setw(10) << connect[4];
          (*output) << std::setw(10) << connect[4];
          (*output) << std::setw(10) << connect[4];
          (*output) << std::setw(10) << connect[4];
          break;

    case Elem::ET_PYRA13:
    case Elem::ET_PYRA14:
          (*output) << std::setw(10) << 20 << std::endl;
          (*output) << std::setw(10) << connect[0];
          (*output) << std::setw(10) << connect[5];

          (*output) << std::setw(10) << connect[1];
          (*output) << std::setw(10) << connect[6];

          (*output) << std::setw(10) << connect[2];
          (*output) << std::setw(10) << connect[7];

          (*output) << std::setw(10) << connect[3];
          (*output) << std::setw(10) << connect[8] << std::endl;

          (*output) << std::setw(10) << connect[9];
          (*output) << std::setw(10) << connect[10];
          (*output) << std::setw(10) << connect[11];
          (*output) << std::setw(10) << connect[12];

          (*output) << std::setw(10) << connect[4];
          (*output) << std::setw(10) << connect[4];

          (*output) << std::setw(10) << connect[4];
          (*output) << std::setw(10) << connect[4] << std::endl;

          (*output) << std::setw(10) << connect[4];
          (*output) << std::setw(10) << connect[4];

          (*output) << std::setw(10) << connect[4];
          (*output) << std::setw(10) << connect[4];
          break;

    case Elem::ET_WEDGE15:
    case Elem::ET_WEDGE18:
          (*output) << std::setw(10) << 15 << std::endl;
          (*output) << std::setw(10) << connect[0];
          (*output) << std::setw(10) << connect[6];

          (*output) << std::setw(10) << connect[1];
          (*output) << std::setw(10) << connect[7];

          (*output) << std::setw(10) << connect[2];
          (*output) << std::setw(10) << connect[8];

          (*output) << std::setw(10) << connect[12];
          (*output) << std::setw(10) << connect[13] << std::endl;
          (*output) << std::setw(10) << connect[14];

          (*output) << std::setw(10) << connect[3];
          (*output) << std::setw(10) << connect[9];

          (*output) << std::setw(10) << connect[4];
          (*output) << std::setw(10) << connect[10];

          (*output) << std::setw(10) << connect[5];
          (*output) << std::setw(10) << connect[11];
          break;

    case Elem::ET_HEXA20:
//          1         9         2        10         3        11         4        12
//         17        18        19        20         5        13         6        14
//          7        15         8        16

          (*output) << std::setw(10) << 20 << std::endl;
          (*output) << std::setw(10) << connect[0];
          (*output) << std::setw(10) << connect[8];

          (*output) << std::setw(10) << connect[1];
          (*output) << std::setw(10) << connect[9];

          (*output) << std::setw(10) << connect[2];
          (*output) << std::setw(10) << connect[10];

          (*output) << std::setw(10) << connect[3];
          (*output) << std::setw(10) << connect[11] << std::endl;

          (*output) << std::setw(10) << connect[16];
          (*output) << std::setw(10) << connect[17];
          (*output) << std::setw(10) << connect[18];
          (*output) << std::setw(10) << connect[19];

          (*output) << std::setw(10) << connect[4];
          (*output) << std::setw(10) << connect[12];

          (*output) << std::setw(10) << connect[5];
          (*output) << std::setw(10) << connect[13] << std::endl;

          (*output) << std::setw(10) << connect[6];
          (*output) << std::setw(10) << connect[14];

          (*output) << std::setw(10) << connect[7];
          (*output) << std::setw(10) << connect[15];
          break;

        default:
          (*output) << connect.GetSize() << std::endl;

          for (UInt ii=0; ii < connect.GetSize(); ii++) { 
            (*output).width(10);
            (*output) << connect[ii];
          }
          break;
        }
        
        (*output) << std::endl;
      } // over elements of group
      elmsgrp++;
    } // over groups
    (*output) << std::setw(6) << -1 << std::endl;
  }

  void SimOutputUnv::ElemType2UnvElemId(Elem::FEType et, UInt & id)
  {
    switch(et)
    {
    case Elem::ET_LINE2:
      id       = 21;
      break;
    case Elem::ET_LINE3:
      id       = 24;
      break;
    case Elem::ET_TRIA3:
      id       = 91;
      break;
    case Elem::ET_TRIA6:
      id       = 92;
      break;
    case Elem::ET_QUAD4:
      id       = 94;
      break;
    case Elem::ET_QUAD8:
    case Elem::ET_QUAD9:
      id       = 95;
      break;
    case Elem::ET_TET4:
      id       = 111;
      break;
    case Elem::ET_TET10:
      id       = 118;
      break;
    case Elem::ET_HEXA8:
      id       = 115;
      break;
    case Elem::ET_HEXA20:
      id       = 116;
      break;
    case Elem::ET_HEXA27:
      id       = 117;
      break;
    case Elem::ET_PYRA5:
      id       = 115;
      break;
    case Elem::ET_PYRA13:
    case Elem::ET_PYRA14:
      id       = 116;
      break;
    case Elem::ET_WEDGE6:
      id       = 112;
      break;
    case Elem::ET_WEDGE15:
    case Elem::ET_WEDGE18:
      id       = 113;
      break;
    default:
      id       = 0;
      break;
    }

  }

  void  SimOutputUnv::Dataset2467()
  {
    //
    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized");

    (*output) << std::setw(6) << -1 << std::endl << std::setw(6)
              << 2467 << std::endl;

    StdVector<Elem*> elems;
    StdVector<UInt> entityNums;

    StdVector<RegionIdType> subdoms;
    ptGrid_->GetVolRegionIds(subdoms);
    StdVector<RegionIdType> surfsubdoms;
    ptGrid_->GetSurfRegionIds(surfsubdoms);

    StdVector<std::string> regionNames;
    ptGrid_->GetRegionNames(regionNames);

    UInt i, j, k,n;
    k = 0;
    for (i=0; i<surfsubdoms.GetSize(); i++) {
      subdoms.Push_back(surfsubdoms[i]);
    }

    StdVector<std::string> names;
    StdVector<UInt> entityTypeCodes;
    StdVector< StdVector<UInt> > entities;

    // Collect region elements
    for (i=0; i<subdoms.GetSize(); i++) {
      UInt numElems = ptGrid_->GetNumElems(subdoms[i]);
      ptGrid_->GetElems(elems,subdoms[i]);

      names.Push_back(regionNames[subdoms[i]]);
      entityTypeCodes.Push_back(8);

      entityNums.Clear();
      for ( j = 0; j < numElems; j++) {  
        entityNums.Push_back(elems[j]->elemNum);
      }

      entities.Push_back(entityNums);      
    }

    StdVector<std::string> nodeNames, elemNames;
    UInt numNamedEntities = 0;

    // Collect named elements
    ptGrid_->GetListElemNames(elemNames);
    numNamedEntities = elemNames.GetSize();

    for( UInt i = 0; i < numNamedEntities; i++ ) {
      
      // get elements
      ptGrid_->GetElemsByName(elems, elemNames[i]);
      UInt numElems = elems.GetSize();

      names.Push_back(elemNames[i]);
      entityTypeCodes.Push_back(8);

      entityNums.Clear();
      for ( j = 0; j < numElems; j++) {  
        entityNums.Push_back(elems[j]->elemNum);
      }

      entities.Push_back(entityNums);      
    }
    
    // Collect named nodes
    ptGrid_->GetListNodeNames(nodeNames);
    numNamedEntities = nodeNames.GetSize();

    for( UInt i = 0; i < numNamedEntities; i++ ) {      
      // get nodes
      entityNums.Clear();
      ptGrid_->GetNodesByName(entityNums, nodeNames[i]);
      
      names.Push_back(nodeNames[i]);
      entityTypeCodes.Push_back(7);
      entities.Push_back(entityNums);      
    }
  
    UInt group=1;
    for (i=0, n=names.GetSize(); i<n; i++) {    
      UInt numEntities = entities[i].GetSize();

      (*output) << std::setw(10) << (group++);
      (*output) << std::setw(10) << 0;
      (*output) << std::setw(10) << 0;
      (*output) << std::setw(10) << 0;
      (*output) << std::setw(10) << 0;
      (*output) << std::setw(10) << 0;
      (*output) << std::setw(10) << 0;
      (*output) << std::setw(10) << numEntities << std::endl;
      (*output) << names[i] << std::endl;

      for ( j = 0; j < numEntities; j++, k++ ) {  
        (*output) << std::setw(10) << entityTypeCodes[i];
        (*output) << std::setw(10) << entities[i][j];
        (*output) << std::setw(10) << 0;
        (*output) << std::setw(10) << 0;
        if( (j % 2) == 1 ) (*output) << std::endl;
      } // over elements of group
      if( (j % 2) == 1 ) (*output) << std::endl;
    } // over groups
    (*output) << std::setw(6) << -1 << std::endl;
  }

  void  SimOutputUnv::NodeElemDataTransient(const UInt dataSetNr,
                                                  const std::string & title, 
                                                  const Vector<Double> & x, 
                                                  const UInt step, 
                                                  const Double time, 
                                                  const UInt nrNodes,
                                                  const UInt nrDofs)
  {
    //
    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized" );

    (*output) << std::setw(6) << -1 << std::endl 
              << std::setw(6) << dataSetNr << std::endl;

    (*output).setf(std::ios::scientific);
    (*output).precision(6);
    (*output).setf(std::ios::uppercase);

    // Standard for scalar values
    UInt dataCharac = 1;
    UInt valsPerNode = 1;
 
    // needed for undocumented value of
    // Dataset 55/56 in record8 , field4
    UInt specDataCharac = 0;

    // Vector type
    if (nrDofs > 1 && nrDofs <= 3) 
      {
        valsPerNode = 3;
        dataCharac = 2;
      } 
    // Tensor
    else if(nrDofs == 6) 
      { 
        valsPerNode = 6;
        specDataCharac = 2;
        dataCharac = 4; // symmetric tensor
        //      dataCharac = 3; // vector with 6 components
        //      dataCharac = 5; // unsymmetric tensor
      }
     
    if(capaOut_) 
    {
      (*output) << " " << title << " step" << std::setw(6) << step <<
        " time   " << time << std::endl;  
      (*output) << std::endl << std::endl << std::endl << std::endl;
    }
    else 
    {
      (*output) << " " << title << std::endl;  
      (*output) << "NONE" << std::endl 
                << "NONE" << std::endl 
                << "NONE" << std::endl 
                << "NONE" << std::endl;
    }
    (*output) << std::setw(10) << 1 << std::setw(10) << 4 << std::setw(10) 
              << dataCharac  << std::setw(10) << specDataCharac
              << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
    (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 
              << std::setw(10) <<
      step << std::endl;
    (*output) << " " << time << std::endl;       

    UInt i,j,n;
    n=nrNodes;;  
    for (i=0; i<n; i++)
      {
     
        (*output) << std::setw(10) << i+1;
        if (dataSetNr == 56)
          (*output) << std::setw(10) << valsPerNode;

        (*output) << std::endl;
     
        // in the universal file either one or three results datas must exist
        if (nrDofs == 2)
          (*output) << 0.0;

        for (j=0; j<nrDofs; j++)
          {
            //std::cerr << "trying to write " << i << ", " << j << std::endl;
            (*output) << std::setw(14) << x[i*nrDofs +j];
          }
     
        (*output) << std::endl;
      }    
    (*output) << std::setw(6) << -1 << std::endl;
  }  

  void SimOutputUnv::NodeElemDataHarmonic(const UInt dataSetNr,
                                                const std::string & title, 
                                                const Vector<Complex> & x, 
                                                const UInt step,
                                                const Double frequency, 
                                                const ComplexFormat format,
                                                const UInt nrNodes,
                                                const UInt nrDofs)
  {
  
    UInt dataCharact = 1;
    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized" );
  
    (*output) << std::setw(6) << -1 << std::endl 
              << std::setw(6) << dataSetNr << std::endl;
  
    (*output).setf(std::ios::scientific);
    (*output).precision(6);
    (*output).setf(std::ios::uppercase);
  
    UInt valsPerNode = 1;
    if (nrDofs > 1)
      {
        dataCharact = 2;
        valsPerNode = 3;
      }
 

    if (format == REAL_IMAG)
      {
        // write out realpart
        if(capaOut_) 
        {
          (*output) << " " << title << " cw realpart" << std::setw(6) <<" frequency   " << frequency << std::endl;  
        } else 
        {        
          (*output) << " " << title << "_RE" << std::endl;  
        }
        
        (*output) << std::endl << std::endl << std::endl << std::endl;
        (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
                  << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
        (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -2 << std::setw(10) <<
          step << std::endl;
        (*output) << " " << frequency << std::endl;       
      
        UInt i,j,n;
        n=nrNodes;
        for (i=0; i<n; i++)
          {
            (*output) << std::setw(10) << i+1 << std::endl;
          
            // in the universal file either one or three results datas must exist
            if (nrDofs == 2)
              (*output) << 0.0;
          
            for (j=0; j<nrDofs; j++)
              {
                //std::cerr << "trying to write " << i << ", " << j << std::endl;
                (*output) << std::setw(14) << x[i*nrDofs +j].real();
              }
          
            (*output) << std::endl;
          }    
        (*output) << std::setw(6) << -1 << std::endl;

        // write out imag part
        (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 55 << std::endl;
        if(capaOut_) 
        {
          (*output) << " " << title << " cw imagpart" << std::setw(6) <<" frequency   " << frequency << std::endl;  
        } else 
        {        
          (*output) << " " << title << "_IM" << std::endl;
        }
        
        (*output) << std::endl << std::endl << std::endl << std::endl;
        (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
                  << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
        (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -12 << std::setw(10) <<
          step << std::endl;
        (*output) << " " << frequency << std::endl;       
      
        for (i=0; i<n; i++)
          {
            (*output) << std::setw(10) << i+1 << std::endl;
          
            // in the universal file either one or three results datas must exist
            if (nrDofs == 2)
              (*output) << 0.0;
          
            for (j=0; j<nrDofs; j++)
              {
                //std::cerr << "trying to write " << i << ", " << j << std::endl;
                (*output) << std::setw(14) << x[i*nrDofs +j].imag();
              }
          
            (*output) << std::endl;
          }    
        (*output) << std::setw(6) << -1 << std::endl;
      }
  
    else if (format == AMPLITUDE_PHASE) {
      // write out amplitude
        if(capaOut_) 
        {
          (*output) << " " << title << " cw amplitude" << std::setw(6) <<" frequency   " << frequency << std::endl;  
        } else {
          (*output) << " " << title << "_AM" << std::endl;
        }        
      (*output) << std::endl << std::endl << std::endl << std::endl;
      (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
                << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
      (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -1 << std::setw(10) <<
        step << std::endl;
      (*output) << " " << frequency << std::endl;       
    
      UInt i,j,n;
      n=nrNodes;
      for (i=0; i<n; i++)
        {
          (*output) << std::setw(10) << i+1 << std::endl;
        
          // in the universal file either one or three results datas must exist
          if (nrDofs == 2)
            (*output) << 0.0;
        
          for (j=0; j<nrDofs; j++)
            {
              //std::cerr << "trying to write " << i << ", " << j << std::endl;
              (*output) << std::setw(14) << std::abs(x[i*nrDofs +j]);
            }
        
          (*output) << std::endl;
        }    
      (*output) << std::setw(6) << -1 << std::endl;
    
      // write out phase
      (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 55 << std::endl;
        if(capaOut_) 
        {
          (*output) << " " << title << " cw phase" << std::setw(6) <<" frequency   " << frequency << std::endl;  
        } else {
          (*output) << " " << title << "_PH" << std::endl;
        }        
      (*output) << std::endl << std::endl << std::endl << std::endl;
      (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
                << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
      (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -11 << std::setw(10) <<
        step << std::endl;
      (*output) << " " << frequency << std::endl;       
    
      for (i=0; i<n; i++)
        {
          (*output) << std::setw(10) << i+1 << std::endl;
        
          // in the universal file either one or three results datas must exist
          if (nrDofs == 2)
            (*output) << 0.0;
        
          for (j=0; j<nrDofs; j++)
            {
              if (abs(x[i*nrDofs +j].imag()) > 1e-16)
                (*output) << std::setw(14) << CPhase(x[i*nrDofs +j]);
              else 
                (*output) << std::setw(14) << 0.0;
            }
        
          (*output) << std::endl;
        }    
      (*output) << std::setw(6) << -1 << std::endl;
    }
  
  }

  void SimOutputUnv::Init(Grid * aptgrid, bool printGridOnly )
  {
    try 
    {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }

    std::string name = fileName_ + "." + formatName_;
    fs::path filePath = dirName_ / name;
    output = NULL;
    output = new std::ofstream(filePath);
    if(!output)
      EXCEPTION("Unv file ' " << name << "' could not be openend!" );

    ptGrid_=aptgrid;
    WriteGrid();

  }

  void SimOutputUnv::RegisterResult( shared_ptr<BaseResult> sol,
                                     UInt saveBegin, UInt saveInc,
                                     UInt saveEnd,
                                     bool isHistory ) {
    
  }

  void SimOutputUnv::BeginStep( UInt stepNum, Double stepVal ) {
    
    actStep_ = stepNum + stepNumOffset_;
    actStepVal_ = stepVal;
    resultMap_.clear();
  }

  void SimOutputUnv::AddResult( shared_ptr<BaseResult>  sol ) {
    resultMap_[sol->GetResultInfo()->resultName].Push_back( sol );
  }
  
  void SimOutputUnv::FinishStep( ) {

    // iterate over all result types
    ResultMapType::iterator it = resultMap_.begin();
    for( ; it != resultMap_.end(); it++ ) {
      
      // check if result is defined on nodes or elements
      ResultInfo & actInfo = *(it->second[0]->GetResultInfo());
      
      const StdVector<shared_ptr<BaseResult> > actResults =
        it->second;

      if(!ValidateNodesAndElements(actInfo)) continue;      
      
      ResultInfo::EntityUnknownType entityType = actInfo.definedOn;
      std::string title;
      std::string ideas5xSolType;
      UInt ideas2414SolType;
      

      SolutionTypeToString( actInfo.resultType, title,
                            ideas5xSolType, ideas2414SolType );
      
      if(!capaOut_) 
      {
        title = ideas5xSolType;
      }
      
      // The ideas2414SolType  is meant to  provide a mapping from  our result
      // types to the I-DEAS types defined for dataset number 2414 which we do
      // not currently write, but may do  so in the future. Here is an example
      // for a scalar and a 6-dof dataset for a three-node triangle:

      //     -1
      //   2414   
      //          3
      //  Sound Pressure 
      //          1
      // NONE  
      // MODEL_SOLUTION_SOLVE
      // LOAD SET 1                   
      // Analysis time was     11-Nov-99       17:09:55   
      // NONE     
      //          1         1         1       301         2         1
      //        -10         0         1         1         1         0         0         0
      //          2         0
      //   0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00
      //   0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00
      //          1
      //   1.00000E-03
      //          2
      //   0.00000E-03
      //          3
      //   2.00000E-03
      //     -1
      //     -1
      //   2414
      //          2
      //  B.C. 1,REACTION FORCE_2,LOAD SET 1
      //          1
      // NONE  
      // MODEL_SOLUTION_SOLVE                                             
      // LOAD SET 1            
      // Analysis time was     11-Nov-99       17:09:55
      // NONE    
      //          1         1         3         9         2         6
      //        -10         0         1         1         1         0         0         0
      //          2         0
      //   0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00
      //   0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00
      //          1
      //   4.43730E-03 -1.75471E-02 -6.28080E-03  0.00000E+00  0.00000E+00  0.00000E+00
      //          2
      //   8.61587E-03 -2.87617E-02 -2.40538E-03  0.00000E+00  0.00000E+00  0.00000E+00
      //          3
      //   8.58183E-03 -2.80713E-02 -9.58135E-04  0.00000E+00  0.00000E+00  0.00000E+00
      //     -1

      // if result could not be mapped, omit it
      if( title == "") {
        WARN("Result '" << actInfo.resultName
             << "' coul not be mappted to unv result type and is "
             << "omitted!");
        continue;
      }
      
      StdVector<std::string> & dofNames = actInfo.dofNames;
      UInt numDofs = dofNames.GetSize();
      
      // Determine type of entity the result is defined on
      UInt mapTo = 0;
      UInt numEntities = 0;
      if( actInfo.definedOn == ResultInfo::NODE ) {
        mapTo = 55;
        numEntities = ptGrid_->GetNumNodes();;
      } else {
        mapTo = 56;
        numEntities = ptGrid_->GetNumVolElems();
      }
      
      if( actResults[0]->GetEntryType() == BaseMatrix::DOUBLE ) {
        Vector<Double> gSol;
        FillGlobalVec<Double>(gSol, actResults, entityType );
        
        // Special case: If result is mechStress, restort entries 
        // according to capa notation
        if( actInfo.resultType == MECH_STRESS ) { 
          SortStresses( gSol, dofNames);
          numDofs = 6;
        }
        
        NodeElemDataTransient( mapTo, title, gSol, actStep_, 
                               actStepVal_ , numEntities ,
                               numDofs );
      } else {
        Vector<Complex> gSol;
        FillGlobalVec<Complex>(gSol, actResults, entityType );

        // Special case: If result is mechStress, restort entries 
        // according to capa notation
        if( actInfo.resultType == MECH_STRESS ) { 
          SortStresses( gSol, dofNames);
          numDofs = 6;
        }
        
        NodeElemDataHarmonic( mapTo, title, gSol, actStep_, 
                              actStepVal_, actInfo.complexFormat,
                              numEntities, numDofs );
      }
    } // over all result types
  }
  
  
  void SimOutputUnv::FinishMultiSequenceStep() {

    // set offset for step value and number to last values
    stepNumOffset_ = actStep_;
  }
  
  
  template<class TYPE>
  void SimOutputUnv::SortStresses( Vector<TYPE> & vec,

                                   const StdVector<std::string>& dofNames ){


    //soring according to capa (unv) notation
    //our notation is Voigt: Txx Tyy Tzz Tyz Txz Txy


    // check size of dofNames:
    // 6: 3D
    // 4: axi
    // 3: plane 

    // ensure, that number of entries is multiples of number of dofs
    assert( vec.GetSize() % dofNames.GetSize() == 0 );
    UInt numDofs = dofNames.GetSize();
    UInt numEntities = (UInt) vec.GetSize() / dofNames.GetSize();

    // CAPA always expects a stress vector with 6 entries
    Vector<TYPE> sorted;
    sorted.Resize( numEntities * 6 );
    
    // a) 3D-case (6 entries for stress vector)
    if( numDofs == 6 ) {
      
      //Txx Txy Tyy Txz Tyz Tzz
      UInt j = 0;
      for( UInt i = 0; i < vec.GetSize(); i+=6 ) {
        sorted[j++] = vec[i];
        sorted[j++] = vec[i+5];
        sorted[j++] = vec[i+1];
        sorted[j++] = vec[i+4];
        sorted[j++] = vec[i+3];
        sorted[j++] = vec[i+2];
      }
    }
  
    // b) axisymmetric-case (4 entries for stress vector)  
    else if( numDofs == 4 ) {
      //Tphiphi 0 Trr 0 Trz Tzz
      
      
      UInt j = 0;
      for( UInt i = 0; i < vec.GetSize(); i+=4 ) {
        sorted[j++] = vec[i+3];
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+0];
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+2];
        sorted[j++] = vec[i+1];
      }
    }

    // c) planestress/strain-case (3 entries for stress vector)  
    else if( numDofs == 3 ) {
      //0 0 Tyy 0 Tyz Tzz
    
      UInt j = 0;
      for( UInt i = 0; i < vec.GetSize(); i+=3 ) {
        sorted[j++] = 0.0;
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+0];
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+2];
        sorted[j++] = vec[i+1];
      } 

    } else {
      EXCEPTION( "Can not resort a stress vector with " << numDofs << " entries" );
    }
   
    // store sorted vector back to original one
    vec = sorted;
  }
  
  void SimOutputUnv::SolutionTypeToString(const SolutionType type,
                                          std::string& capaSolType,
                                          std::string& ideas5xSolType,
                                          UInt& ideas2414SolType) const
  {
    capaSolType = "";
    ideas5xSolType = "UNK";
    ideas2414SolType = 0;
    switch (type)
      {
        // MECHANICS
      case MECH_DISPLACEMENT:
        capaSolType = "displacement";
        ideas5xSolType = "U";
        ideas2414SolType = 8;
        break;
      case MECH_ACCELERATION:
        capaSolType = "acceleration";
        ideas5xSolType = "A";
        ideas2414SolType = 12;
        break;
      case MECH_VELOCITY:
        capaSolType = "velocity";
        ideas5xSolType = "V";
        ideas2414SolType = 11;
        break;
      case MECH_FORCE:
        ideas5xSolType = "F";
        ideas2414SolType = 9;
        break;
      case MECH_STRESS:
        capaSolType = "stress";
        ideas5xSolType = "STRESS";
        ideas2414SolType = 2;
        break;
      case MECH_STRAIN:
        EXCEPTION("Not implemented" );
        ideas5xSolType = "STRAIN";
        ideas2414SolType = 3;
        break;
      case MECH_RHS_LOAD:
        capaSolType = "mechRhsLoad";
        ideas5xSolType = "MECHRHS";
        ideas2414SolType = 4;
        break;

        // ELECTRIC
      case ELEC_POTENTIAL:
        capaSolType = "electric potential";
        ideas5xSolType = "VOLT";
        ideas2414SolType = 94;
        break;
      case ELEC_FIELD_INTENSITY:
        capaSolType = "electric field";
        ideas5xSolType = "EFIELD";
        ideas2414SolType = 95;
        break;
      case ELEC_FORCE_VWP: 
        EXCEPTION("Not implemented" );
        break;
      case ELEC_CHARGE:
        capaSolType = "electric charge";
        ideas5xSolType = "ECHARGE";
        break;
      case ELEC_FLUX_DENSITY:
        EXCEPTION("Not implemented");
        break; 
      case ELEC_ENERGY:
        EXCEPTION("Not implemented");
        break;
      case ELEC_RHS_LOAD:
        capaSolType = "elecRhsLoad";
        ideas5xSolType = "ERHS";
        break;
      case SMOOTH_DISPLACEMENT:
        capaSolType = "displacement";
        ideas5xSolType = "USMOOTH";
        break;

        // ACOUSTIC
      case ACOU_POTENTIAL:
        capaSolType = "fluid potential";
        ideas5xSolType = "APOT";
        ideas2414SolType = 301;
        break;
      case ACOU_RHS_LOAD:
        capaSolType = "fluid potential";
        ideas5xSolType = "ALHRHS";
        ideas2414SolType = 94;
        break;
      case ACOU_PRESSURE:
        //       warnMsg = "Due to the restrictions in the .unv file format, the ";
        //       warnMsg += "acoustic pressure is written as acoustic (fluid) potential!";
        //       WARN(warnMsg.c_str());
        capaSolType = "fluid potential";
        ideas5xSolType = "APRES";
        ideas2414SolType = 301;
        break;
      case ACOU_FORCE:
        EXCEPTION("Not implemented" );
        break;
      case ACOU_POTENTIAL_DERIV_1:
        capaSolType = "fluid potential, 1st deriv.";
        ideas5xSolType = "APOTD1";
        ideas2414SolType = 94;
        break;
      case ACOU_POTENTIAL_DERIV_2:
        capaSolType = "fluid potential, 2nd deriv.";
        ideas5xSolType = "APOTD2";
        break;

        // MAGNETIC
      case MAG_POTENTIAL:
        capaSolType = "mag. vector potential";
        ideas5xSolType = "MAGPOT";
        break;
      case MAG_FLUX_DENSITY:
        capaSolType = "mag. flux density";
        ideas5xSolType = "MAGFLUX";
        break;
      case MAG_EDDY_CURRENT_DENSITY:
        capaSolType = "eddy current";
        ideas5xSolType = "MAGEDDY";
        break;
      case MAG_FORCE_VWP:
        EXCEPTION("Not implemented");
        break;
      case MAG_FORCE_LORENTZ_DENSITY:
        EXCEPTION("Not implemented");
        break;
      case MAG_ENERGY:
        EXCEPTION("Not implemented");
        break;
      case MAG_RHS_LOAD:
        capaSolType = "magRhsLoad";
        ideas5xSolType = "MAGRHS";
        break;

        // HEAT
      case HEAT_TEMPERATURE:
        capaSolType = "temperature";
        ideas5xSolType = "TEMP";
        ideas2414SolType = 5;
        break;
      case HEAT_RHS_LOAD:
        capaSolType = "temperatureRhsLoad";
        ideas5xSolType = "TEMPRHS";
        break;

        // FLUID
      case FLUIDMECH_VELOCITY:
        capaSolType = "velocity";
        ideas5xSolType = "VEL";
        ideas2414SolType = 11;
        break;
      case FLUIDMECH_PRESSURE:
        capaSolType = "pressure";
        ideas5xSolType = "PRES";
        ideas2414SolType = 117;
        break;
      case MEAN_FLUIDMECH_VELOCITY:
        capaSolType = "mean_veclocity";
        ideas5xSolType = "MEANVEL";
        ideas2414SolType = 53;
        break;

      default:
        capaSolType = "";
      }
  }




  // void  SimOutputUnv::WriteDataOnCell(const Vector<Double> & sol, const UInt step, const Double time, const std::string title, const Integer nrDofs)
  // {
  //   std::string aux;
  //   if (title == "elecfield") {
  //     aux="electric field";
  //     Dataset56(aux, sol, step+1, time);
  //   }
  //   else 
  //     Dataset56(title, sol, step+1, time, nrDofs);
  //     //WARN("This cell-data is not printed, since this type is not supported by Capapost");
  // }

}
