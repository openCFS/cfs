#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acoustic3dPDE.hh"
#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>


namespace CoupledField
{

Acoustic3dPDE::Acoustic3dPDE(Grid * aptgrid, BCs *aptbcs, Material *ptMaterial, TimeFunc *aptTimeFunc, 
		     FileType *aptFileType, WriteResults *aptOut)
:AcousticPDE(aptgrid,aptbcs,ptMaterial,aptTimeFunc,aptFileType,aptOut)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::Acoustic3dPDE " << std::endl;
#endif

  pdename_    ="Acoustic3d";

  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);

}


void Acoustic3dPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::SetupMatrices" << std::endl;
#endif

  Matrix<Double> elemmat;
  Matrix<Double> elemmatbnd; // This is a smaller matrix since it is just for linear 1D elements.
  Point<3> * ptCoord;

  BaseElem * ptEl;

  Vector<Double> coeffm, coeffst, coeffdamp;
  CalcCoeff(coeffm, coeffst, coeffdamp);

  Vector<Integer> connecth;
  std::vector<Elem*> elemssd;

  Integer i, j;
  for (i=0; i<subdoms_.size(); i++)
    {
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j< elemssd.size(); j++)
	{
	  ptEl=elemssd[j]->ptElem;

	  BaseForm<3> * bilinear_stiff = new LaplaceInt<3>(ptEl,1);
	  BaseForm<3> * bilinear_mass  = new MassInt<3>(ptEl,1);

	  Integer ii;
	  Integer elsize=(elemssd[j]->connect).size();
	  connecth.Resize(elsize);
	  for (ii=0; ii<elsize; ii++)
	    connecth[ii]=(elemssd[j]->connect)[ii];

	  ptCoord=new Point<3>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);
	  elemmat*=coeffst[i];

#ifdef DEBUG
	  (*debug) << "Connection array  " << std::endl;
	  (*debug)  << connecth  << std::endl;
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   j << std::endl;
	  (*debug) << elemmat << std::endl;
#endif

	  algsys_->SetElementMatrix(elemmat.getinarray(), connecth.get(), connecth.size(), STIFFNESS);

	  // mass part
	  bilinear_mass->CalcElemMatrix(ptCoord, elemmat);
	  elemmat*=coeffm[i];

#ifdef DEBUG
	  (*debug) << "Massmatrix, ElementNumber  " << j << std::endl;
 	  (*debug) << elemmat << std::endl;
#endif

	  algsys_->SetElementMatrix(elemmat.getinarray(), connecth.get(), connecth.size(), MASS);

	  delete bilinear_stiff;
	  delete bilinear_mass;
	  delete [] ptCoord;
	}
    }

  // BEGIN DAMPING MATRIX PART
  if (with_absBCs_) {
  std::vector<Elem*>  DomainBnd;
  conf->getliststr("bnd_for_absBCs",bnd_absBCs_,"Acoustic");
  DomainBnd=ptBCs_->getFacesBC(bnd_absBCs_[0],level);

  for (j=0; j< DomainBnd.size(); j++)

    {
      ptEl=DomainBnd[j]->ptElem;

      BaseForm<3> * linear_damp = new DampInt<3>(ptEl,1);

      Integer ii;
      Integer elsize=(DomainBnd[j]->connect).size();
      connecth.Resize(elsize);
      for (ii=0; ii<elsize; ii++)
	connecth[ii]=(DomainBnd[j]->connect)[ii];

      ptCoord=new Point<3>[connecth.size()];
      ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

      linear_damp->CalcElemMatrix(ptCoord, elemmatbnd);
      elemmatbnd*=coeffdamp[0];

#ifdef DEBUG
      (*debug) << "Connection array  " << std::endl;
      (*debug)  << connecth  << std::endl;
      (*debug) << "Dampingmatrix, ElementNumber  " <<   j << std::endl;
      (*debug) << elemmatbnd << std::endl;
#endif

      algsys_->SetElementMatrix(elemmatbnd.getinarray(), connecth.get(), connecth.size(), DAMPING);

      delete linear_damp;
      delete [] ptCoord;
    }
}      // END DAMPING MATRIX PART


#ifdef TRACE
  (*trace) << "Leaving Acoustic3dPDE::SetupMatrices" << std::endl;
#endif
}

Acoustic3dPDE::~Acoustic3dPDE()
{
  ;
}



} // end of namespace

