#include <vector>
#include <string>
#include <fstream>

#include "Elements/elements_header.hh"
#include "grid.hh"
#include "DataInOut/conffile.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{

Grid::Grid(FileType * aptFileType)
{
#ifdef TRACE
  (*trace) << " entering Grid::Grid " << std::endl;
#endif

  ptFileType = aptFileType;

  ptQ   = new Quad1FE();
  ptQ2   = new Quad2FE();
  ptTet = new Tetra1FE();
  ptL1  = new Line1FE();
  ptL2  = new Line2FE();
  ptTr1 = new Triangle1FE();
  ptTr2 = new Triangle2FE();
  ptHexa = new Hexa1FE();
  ptPyra = new Pyra1FE();
  ptWedge = new Wedge1FE();

  lastlevel_=0;

#ifndef XMLPARAMS
  conf->getsubdom(listSD_);
#else
  params->GetList( "name", listSD_, "domain", "region" );
#endif

}



  
void Grid::SetIntTypeAllElems(IntegrationType aIntType)
{
#ifdef TRACE
  (*trace) << " entering Grid::SetIntTypeAllElems " << std::endl;
#endif

  ptQ   -> SetIntegrationType(aIntType);
  ptQ2   -> SetIntegrationType(aIntType);
  ptTet -> SetIntegrationType(aIntType);
  ptL1  -> SetIntegrationType(aIntType);
  ptL2  -> SetIntegrationType(aIntType);
  ptTr1 -> SetIntegrationType(aIntType);
  ptTr2 -> SetIntegrationType(aIntType);
  ptHexa-> SetIntegrationType(aIntType);
  ptPyra-> SetIntegrationType(aIntType);
  ptWedge-> SetIntegrationType(aIntType);

  ptQ   ->Init();
  ptQ2   ->Init();
  ptTet ->Init();
  ptL1  ->Init();
  ptL2  ->Init();
  ptTr1 -> Init();
  ptTr2 -> Init();
  ptHexa-> Init();
  ptPyra-> Init();
  ptWedge-> Init();

}




Grid::~Grid()
{

 if (ptQ)   delete ptQ;
 if (ptQ2)   delete ptQ2;
 if (ptTet) delete ptTet;
 if (ptL1)  delete ptL1;
 if (ptL2)  delete ptL2;
 if (ptTr1)  delete ptTr1;
 if (ptTr2)  delete ptTr2;
 if (ptHexa) delete ptHexa;
 if (ptPyra) delete ptPyra;
 if (ptWedge) delete ptWedge;

}

}
