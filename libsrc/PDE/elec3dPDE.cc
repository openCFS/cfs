#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "elec3dPDE.hh"
#include <Domain/bcs.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>
 
namespace CoupledField
{

Elec3dPDE::Elec3dPDE(Grid *aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, 
		     FileType *aptFileType, WriteResults *aptOut)
:ElecPDE(aptgrid, aptbcs, aptTimeFunc, aptFileType, aptOut)
{
#ifdef TRACE
  (*trace) << "entering Elec3dPDE::Elec3dPDE " << std::endl;
#endif

  pdename_    = "electric3d";
  pdematerialclass_ = "piezo";

  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);
}


void Elec3dPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Elec3dPDE::SetupMatrices" << std::endl;
#endif
  
  Matrix<Double> elemmat;
  Matrix<Double> ptCoord;

  BaseFE * ptElem;


 //reads eps33 (matrix notation starts with 0)
  Double eps33 = materialData_->GetPermittivity(2,2);

  Vector<Integer> connecth;

  Integer i, j;
  for (i=0; i<subdoms_.size(); i++)
    {
      std::vector<Elem*> elemssd;
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j < elemssd.size(); j++)
	{  
	  ptElem=elemssd[j]->ptElem;

	  BaseForm * bilinear_stiff = new LaplaceInt(ptElem, eps33);

	  connecth=elemssd[j]->connect;

	  ptgrid_->GetCoordNodesElemMat(connecth, ptCoord, level);

	  // stiffness part
	  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	  Integer jj;
	  (*debug) << " Pt coords of element " << j << std::endl;
	  (*debug) << ptCoord <<  std::endl;
	  (*debug) << "Stiffnessmatrix, ElementNumber  " << j << std::endl;
	  (*debug) << elemmat << std::endl;
#endif

	  algsys_->SetElementMatrix(elemmat.getinarray(), connecth.get(), connecth.size(), SYSTEM);

	  delete bilinear_stiff;
	}     
    }

#ifdef TRACE
  (*trace) << "Leaving Elec3dPDE::SetupMatrices" << std::endl;
#endif
}

Elec3dPDE::~Elec3dPDE()
{
#ifdef TRACE
  (*trace) << "entering Elec3dPDE::~Elec3dPDE " << std::endl;
#endif
}

} // end of namespace


