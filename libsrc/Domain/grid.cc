#include <vector>
#include <string>
#include <fstream>

#include "elements_header.hh"
#include "line.hh"
#include "grid.hh"
#include "conffile.hh"

namespace CoupledField
{

Grid::Grid(FileType * aptFileType)
{
#ifdef TRACE
  (*trace) << " entering Grid::Grid " << std::endl;
#endif

  ptFileType=aptFileType;
  ptQ=new Quad1();
  ptTr=new Triangle1();
  ptTet=new Tetrahedral1();
  ptL1=new Line();
  ptHexa=new Hexahedral1();
  
  lastlevel_=0;

  conf->getsubdom(listSD_);
}

Grid::~Grid()
{
 if (ptQ) delete ptQ;
 if (ptTr) delete ptTr;
 if (ptTet) delete ptTet;
 if (ptL1) delete ptL1;
 if (ptHexa) delete ptHexa;
}

}
