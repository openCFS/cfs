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

Elec2dPDE::Elec2dPDE(Grid *aptgrid, BCs *aptbcs, Material *ptMaterial, TimeFunc *aptTimeFunc, 
		     FileType *aptFileType, WriteResults *aptOut)
:ElecPDE(aptgrid, aptbcs, ptMaterial, aptTimeFunc, aptFileType, aptOut)

{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::Elec2dPDE " << std::endl;
#endif

  pdename_    = "Electric2d";

  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);
}


void Elec2dPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::SetupMatrices" << std::endl;
#endif
  
  Matrix<Double> elemmat;
  Point<2> * ptCoord;

  BaseElem * ptElem;

  if (InfoPrint)
    (*infofile) << " ------------------------- Element matrices --------------- " << std::endl;

  Vector<Double> coeffst;
  
  CalcCoeff(coeffst);  

  Integer i, j;
  for (i=0; i<subdoms_.size(); i++)
    {
      std::vector<Elem*> elemssd;
   
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j < elemssd.size(); j++)
	{  
	  Vector<Integer> connecth;
 
	  ptElem=elemssd[j]->ptElem;

	  BaseForm<2> * bilinear_stiff = new LaplaceInt<2>(ptElem,1);

	  connecth=elemssd[j]->connect;

	  ptCoord=new Point<2>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);
	  elemmat*=coeffst[i];
	  
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
