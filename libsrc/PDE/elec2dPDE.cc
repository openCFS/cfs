#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "elec2dPDE.hh"
#include <DataInOut/Unverg/outUnverg.hh>
#include <Forms/forms_header.hh>
#include <AlgebraicSystem/LinAlg/linsystem.hh>
#include <Estimator/spaceerror.hh>
#include <multigrid.hh>

namespace CoupledField
{

Elec2dPDE::Elec2dPDE(Grid *aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		     WriteResults *aptOut)
:ElecPDE(aptgrid, aptbcs, aptTimeFunc, aptFileType, aptOut)

{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::Elec2dPDE " << std::endl;
#endif

  pdename_    = "electric2d";
  pdematerialclass_ = "piezo";

  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);
}


void Elec2dPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::SetupMatrices" << std::endl;
#endif
  
  Matrix<Double> elemmat;  
  Matrix<Double> ptCoord;

  BaseFE * ptElem;

  if (InfoPrint)
    (*infofile) << " ------------------------- Element matrices --------------- " << std::endl;

  Vector<Integer> connecth;  

  //reads eps33 (matrix notation starts with 0)
  Double eps0 = materialData_->GetPermittivity(2,2);

  Integer i, j;


  for (i=0; i<subdoms_.size(); i++)
    {
      std::vector<Elem*> elemssd;
   
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j < elemssd.size(); j++)
	{  
	  ptElem=elemssd[j]->ptElem;

	  BaseForm * bilinear_stiff = new LaplaceInt(ptElem, materialData_[i]);

	  connecth=elemssd[j]->connect;
	  
	  ptgrid_->GetCoordNodesElemMat(connecth, ptCoord, level);

	  // stiffness part
	  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);
	  elemmat *= eps0;
	  
#ifdef DEBUG
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif

	  algsys_->SetElementMatrix(elemmat.getinarray(), connecth.get(), connecth.size(), SYSTEM);
	  
	  delete bilinear_stiff;
	  

	}  
      
    }

#ifdef TRACE
  (*trace) << "Leaving Elec2dPDE::SetupMatrices" << std::endl;
#endif
}


Elec2dPDE::~Elec2dPDE()
{
#ifdef TRACE
  (*trace) << " entering Elec2dPDE::~Elec2dPDE " << std::endl;
#endif

}

} // end of namespace
