#include <fstream>
#include <iostream>
#include <string>

#include "basedriver.hh"

#include "DataInOut/AnsysOut/outAnsys.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "DataInOut/Unverg/outUnverg.hh"
#include <DataInOut/writeresults.hh>
namespace CoupledField
{

BaseDriver :: BaseDriver(Domain * adomain)
{
  ENTER_FCN( "BaseDriver::BaseDriver" );

  ptdomain_ = adomain;

  nummeshes_=0;
}

BaseDriver :: ~BaseDriver()
{
  ENTER_FCN( "BaseDriver::~BaseDriver" );
  delete ptdomain_;
}

// for computation with adaptivity
Boolean BaseDriver::printMeshesOrNot()
{
  ENTER_FCN( "BaseDriver::DefinePrintMeshesOrNot" );
  
  Error("Currently not working, need change to XML-Standard");
  
  Boolean meshesInfo=FALSE;
//   std::string typeForMeshesInfo;
//   if (conf->ifget("print_meshes",typeForMeshesInfo))
//     if (typeForMeshesInfo=="ansys")
//       {
// 	system("rm -rf ./mesh_ansys");
// 	//	ptMeshes_=new WriteResultsAnsys("mesh");
// 	meshesInfo = TRUE;
//       }
//    else
//      if (typeForMeshesInfo=="gmv")
//        {
// 	 system("rm -rf ./mesh_gmv");  
// 	 ptMeshes_=new WriteResultsGMV("mesh");
// 	 meshesInfo = TRUE;
//        }
//      else
//        if (typeForMeshesInfo=="unverg")
//        {
// 	 //	 system("rm -rf ./mesh_unv");  
// 	 ptMeshes_=new WriteResultsUnverg("mesh");
// 	 meshesInfo = TRUE;
//        }
  
  return meshesInfo;
}

void BaseDriver :: PrintSeqMeshes()
{
  ENTER_FCN( "BaseDriver::PrintSeqMeshes" );

   if (nummeshes_) {
     ptMeshes_->OpenFile(nummeshes_);
   }

   Integer level=0;
   ptMeshes_->WriteGrid(level);

   nummeshes_++;
}

}
