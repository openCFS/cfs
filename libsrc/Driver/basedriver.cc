#include <fstream>
#include <iostream>
#include <string>

#include "basedriver.hh"

#include "DataInOut/AnsysOut/outAnsys.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/conffile.hh"
#include <DataInOut/writeresults.hh>
namespace CoupledField
{

BaseDriver :: BaseDriver(Domain * adomain)
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::BaseDriver" << std::endl;
#endif
  ptdomain_ = adomain;

  nummeshes_=0;
 
}

BaseDriver :: ~BaseDriver()
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::~BaseDriver" << std::endl;
#endif

}

// for computation with adaptivity
Boolean BaseDriver::printMeshesOrNot()
{
#ifdef TRACE
  (*trace) << " entering BaseDriver :: DefinePrintMeshesOrNot \n";
#endif
  
  Boolean meshesInfo=FALSE;
  std::string typeForMeshesInfo;
  if (conf->ifget("print_meshes",typeForMeshesInfo))
    if (typeForMeshesInfo=="ansys")
      {
	system("rm -rf ./mesh_ansys");
	//	ptMeshes_=new WriteResultsAnsys("mesh");
	meshesInfo = TRUE;
      }
   else
     if (typeForMeshesInfo=="gmv")
       {
	 system("rm -rf ./mesh_gmv");  
	 ptMeshes_=new WriteResultsGMV("mesh");
	 meshesInfo = TRUE;
       }
     else
       if (typeForMeshesInfo=="unverg")
       {
	 //	 system("rm -rf ./mesh_unv");  
	 ptMeshes_=new WriteResultsUnverg("mesh");
	 meshesInfo = TRUE;
       }
  
  return meshesInfo;
}

void BaseDriver :: PrintSeqMeshes()
{
#ifdef TRACE
  (*trace) << "entering BaseDriver :: PrintSeqMeshes" << std::endl;
#endif

   if (nummeshes_) {
     ptMeshes_->OpenFile(nummeshes_);
   }

   Integer level=0;
   ptMeshes_->WriteGrid(level);

   nummeshes_++;
}

}
