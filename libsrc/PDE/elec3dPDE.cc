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

Elec3dPDE::Elec3dPDE(Grid *aptgrid, BCs *aptbcs, Material *ptMaterial, TimeFunc *aptTimeFunc, 
		     FileType *aptFileType, WriteResults *aptOut)
:ElecPDE(aptgrid, aptbcs, ptMaterial, aptTimeFunc, aptFileType, aptOut)
{
#ifdef TRACE
  (*trace) << "entering Elec3dPDE::Elec3dPDE " << std::endl;
#endif

  pdename_    = "Electric3d";
  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);
}


void Elec3dPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Elec3dPDE::SetupMatrices" << std::endl;
#endif
  
  Matrix<Double> elemmat;
  Point<3> * ptCoord;

  BaseElem * ptElem;

  Vector<Double> coeffst;
  CalcCoeff(coeffst);  

  Vector<Integer> connecth;
  std::vector<Elem*> elemssd;

  Integer i, j;
  for (i=0; i<subdoms_.size(); i++)
    {
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j < elemssd.size(); j++)
	{  
	  ptElem=elemssd[j]->ptElem;

	  BaseForm<3> * bilinear_stiff = new LaplaceInt<3>(ptElem,1);

	  connecth=elemssd[j]->connect;

	  ptCoord=new Point<3>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);
	  elemmat*=coeffst[i];

	  if (InfoPrint)
	    (*infofile) << elemmat << std::endl;

#ifdef DEBUG
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif

	  algsys_->SetElementMatrix(elemmat.getinarray(), connecth.get(), connecth.size(), SYSTEM);

	  delete bilinear_stiff;
	  delete [] ptCoord;
	}     
    }

#ifdef TRACE
  (*trace) << "Leaving Elec3dPDE::SetupMatrices" << std::endl;
#endif
}

Elec3dPDE::~Elec3dPDE()
{
 ;
}

} // end of namespace


