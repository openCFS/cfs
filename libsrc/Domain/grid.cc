#include <string>
#include <fstream>

#include "Elements/elements_header.hh"
#include "DataInOut/filetype.hh"
#include "Domain/elem.hh"
#include "grid.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{

  Grid::Grid(FileType * aptFileType)
  {
    ENTER_FCN( "Grid::Grid" );

    ptFileType = aptFileType;

    ptQ1     = new Quad1FE();
    ptQ2     = new Quad2FE();
    ptTet1   = new Tetra1FE();
    ptL1     = new Line1FE();
    ptL2     = new Line2FE();
    ptTr1    = new Triangle1FE();
    ptTr2    = new Triangle2FE();
    ptHexa1  = new Hexa1FE();
    ptHexa2  = new Hexa2FE();
    ptPyra1  = new Pyra1FE();
    ptWedge1 = new Wedge1FE();
    ptWedge2 = new Wedge2FE();


  }



  
  void Grid::SetIntTypeAllElems(IntegrationType aIntType)
  {
    ENTER_FCN( "Grid::SetIntTypeAllElems" );

    ptQ1    -> SetIntegrationType(aIntType);
    ptQ2    -> SetIntegrationType(aIntType);
    ptTet1  -> SetIntegrationType(aIntType);
    ptL1    -> SetIntegrationType(aIntType);
    ptL2    -> SetIntegrationType(aIntType);
    ptTr1   -> SetIntegrationType(aIntType);
    ptTr2   -> SetIntegrationType(aIntType);
    ptHexa1 -> SetIntegrationType(aIntType);
    ptHexa2 -> SetIntegrationType(aIntType);
    ptPyra1 -> SetIntegrationType(aIntType);
    ptWedge1-> SetIntegrationType(aIntType);
    ptWedge2-> SetIntegrationType(aIntType);

    ptQ1    ->Init();
    ptQ2    ->Init();
    ptTet1  ->Init();
    ptL1    ->Init();
    ptL2    ->Init();
    ptTr1   -> Init();
    ptTr2   -> Init();
    ptHexa1 -> Init();
    ptHexa2 -> Init();
    ptPyra1 -> Init();
    ptWedge1-> Init();
    ptWedge2-> Init();

  }




  Grid::~Grid()
  {
    ENTER_FCN( "Grid::~Grid" );
    if (ptQ1)     delete ptQ1;
    if (ptQ2)     delete ptQ2;
    if (ptTet1)   delete ptTet1;
    if (ptL1)     delete ptL1;
    if (ptL2)     delete ptL2;
    if (ptTr1)    delete ptTr1;
    if (ptTr2)    delete ptTr2;
    if (ptHexa1)  delete ptHexa1;
    if (ptHexa2)  delete ptHexa2;
    if (ptPyra1)  delete ptPyra1;
    if (ptWedge1) delete ptWedge1;
    if (ptWedge2) delete ptWedge2;
  
  }




  RegionIdType Grid::RegionNameToId( const std::string & regionName ) {
    ENTER_FCN( "Grid::RegionNameToId" );

    RegionIdType ret = NO_REGION_ID;

    if (regionName == "all" ) {
      ret= ALL_REGIONS;
    } else {
      ret =  regionNames_.Find(regionName);
      if (ret == -1 ) {
        (*error) << "The region with name '" << regionName 
                 << "' is not contained in the grid!";
        Error( __FILE__, __LINE__ );
      }
    }
    return ret;
  }

  void Grid::RegionIdToName( StdVector<std::string> & regionNames,
                             const StdVector<RegionIdType> & regionId ) {
    ENTER_FCN( "Grid::RegionIdToName" );
  
    regionNames.Resize( regionId.GetSize() );
    for (UInt i=0; i<regionId.GetSize(); i++ ) {
      if ( regionId[i] == ALL_REGIONS )
        regionNames[i] = "all";
      else
        regionNames[i] = regionNames_[regionId[i]];
    }
  }

  void Grid::RegionNameToId( StdVector<RegionIdType> & regionIds,
                             const StdVector<std::string> 
                             & regionNames ) {
    ENTER_FCN( "Grid::RegionNameToId" );

    RegionIdType ret = NO_REGION_ID;
    
    regionIds.Resize( regionNames.GetSize() );
    for (UInt i=0; i<regionNames.GetSize(); i++ ) {
    
      if (regionNames[i] == "all" ) {
        ret = ALL_REGIONS;
      } else {
        ret = regionNames_.Find(regionNames[i]);
        if ( ret == -1 ) {
          (*error) << "The region with name '" << regionNames_[i]
                   << "' is not contained in the grid!";
          Error( __FILE__, __LINE__ );
        }
      }
      regionIds[i] = ret;
    }
  }
  
  std::string Grid::RegionIdToName( const RegionIdType regionId ) {
    ENTER_FCN( "Grid::RegionIdToName" );
  
    if ( regionId == ALL_REGIONS )
      return "all";
    else
      return regionNames_[regionId];
  }

} // end of namespace
