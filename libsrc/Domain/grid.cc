#include <vector>
#include <string>
#include <fstream>

#include "Elements/elements_header.hh"
//#include "Elements/1D/line.hh"
#include "grid.hh"
#include "DataInOut/conffile.hh"

namespace CoupledField
{

Grid::Grid(FileType * aptFileType)
{
#ifdef TRACE
  (*trace) << " entering Grid::Grid " << std::endl;
#endif

  ptFileType = aptFileType;
  ptQ   = new Quad1FE();
  ptTet = new Tetra1FE();

//   ptTr=new Triangle1();
//   ptL1=new Line();
//   ptHexa=new Hexahedral1();
  
  lastlevel_=0;

  conf->getsubdom(listSD_);
}

Grid::~Grid()
{
 if (ptQ)   delete ptQ;
 if (ptTet) delete ptTet;

//  if (ptTr) delete ptTr;
//  if (ptL1) delete ptL1;
//  if (ptHexa) delete ptHexa;
}

}
