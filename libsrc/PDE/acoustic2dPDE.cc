#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acoustic2dPDE.hh" 
#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>
#include <Forms/massInt.hh>

namespace CoupledField
{

  Acoustic2dPDE::Acoustic2dPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
			       WriteResults *aptOut)
    :AcousticPDE(aptgrid,aptbcs,aptTimeFunc,aptFileType,aptOut)
  {
#ifdef TRACE
    (*trace) << "entering Acoustic2dPDE::Acoustic2dPDE " << std::endl;
#endif

    pdename_    ="acoustic2d";
    pdematerialclass_ = "fluid";

    conf->getsubdompde(subdoms_,pdename_);
    ReadBCs(pdename_);

    preComputeRHS();
  }


  void Acoustic2dPDE::SetupMatrices(const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering Acoustic2dPDE::SetupMatrices" << std::endl;
#endif

    Matrix<Double> elemmat;

    // This is a smaller matrix since it is just for linear 1D elements.
    Matrix<Double> elemmatbnd;

    BaseFE * ptEl;

    //waiting for material class to be ready
    Double coeffstiff = 1e3;
    Double coeffmass  = 1e3*1e3/2.25e9;

  //  coeffmass[i]  = density*density/compress;
  //  coeffstiff[i] = density;
  //  Double coeffdamp  = density/((sqrt(compress/density)));

    Vector<Integer> connecth;
    std::vector<Elem*> elemssd;

    Integer i, j;
    for (i=0; i<subdoms_.size(); i++)
      {
	ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

	for (j=0; j< elemssd.size(); j++)
	  {
	    ptEl = elemssd[j]->ptElem;

	    BaseForm * bilinear_mass  = new MassInt(ptEl, materialData_[i]);
	    BaseForm * bilinear_stiff = new LaplaceInt(ptEl, materialData_[i]);

	    connecth=elemssd[j]->connect;

	    Matrix<Double> ptCoord;
	    ptgrid_->GetCoordNodesElemMat(connecth, ptCoord, level); 

	    // stiffness part
	    bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);
	    elemmat *= coeffstiff;

#ifdef DEBUG
	    (*debug) << "Connection array  " << std::endl;
	    (*debug)  << connecth  << std::endl;
	    (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	    (*debug) << elemmat << std::endl;
#endif     

	    algsys_->SetElementMatrix(elemmat.getinarray(), connecth.get(), connecth.size(), STIFFNESS);

	    // mass part
	    bilinear_mass->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	    (*debug) << "Massmatrix, ElementNumber  " << i << std::endl;
	    (*debug) << elemmat << std::endl;
#endif
      
	    algsys_->SetElementMatrix(elemmat.getinarray(), connecth.get(), connecth.size(), MASS);
  
	    delete bilinear_stiff;
	    delete bilinear_mass;
	  }
      }

    // BEGIN DAMPING MATRIX PART
    //   if (with_absBCs_) 
    //     {
    //       std::vector<Elem*>  DomainBnd;
    //       conf->getliststr("bnd_for_absBCs",bnd_absBCs_,pdename_);
    //       DomainBnd=ptBCs_->getEdgesBC(bnd_absBCs_[0],level);

    //       for (j=0; j< DomainBnd.size(); j++)
    // 	{
    // 	  ptEl=DomainBnd[j]->ptElem;

    // 	  BaseForm<2> * linear_damp = new DampInt<2>(ptEl,1);

    // 	  Integer ii;
    // 	  Integer elsize=(DomainBnd[j]->connect).size();
    // 	  connecth.Resize(elsize);
    // 	  for (ii=0; ii<elsize; ii++)
    // 	    connecth[ii]=(DomainBnd[j]->connect)[ii];

    // 	  ptCoord=new Point<2>[connecth.size()];
    // 	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

    // 	  linear_damp->CalcElemMatrix(ptCoord, elemmatbnd);
    // 	  elemmatbnd*=coeffdamp[0];

    // #ifdef DEBUG
    // 	  (*debug) << "Connection array  " << std::endl;
    // 	  (*debug)  << connecth  << std::endl;
    // 	  (*debug) << "Dampingmatrix, ElementNumber  " <<   j << std::endl;
    // 	  (*debug) << elemmatbnd << std::endl;
    // #endif

    // 	  algsys_->SetElementMatrix(elemmatbnd.getinarray(), connecth.get(), connecth.size(), DAMPING);

    // 	  delete linear_damp;
    // 	  delete [] ptCoord;
    // 	}
    //    }
    // END DAMPING MATRIX PART

#ifdef TRACE
    (*trace) << "Leaving Acoustic2dPDE::SetupMatrices" << std::endl;
#endif
  }


  void Acoustic2dPDE::preComputeRHS()
  {
#ifdef TRACE
    (*trace) << "entering Acoustic2dPDE::preComputeRHS" << std::endl;
#endif

    SetRHSFnc = FALSE;

    std::string isthererhs="no";
    conf->ifget("load_force",isthererhs,pdename_);

    if (isthererhs=="yes" ) {
      conf->getliststr("rhs_surfaces",rhs_surfaces_,pdename_);
      std::string rhs;
      conf->get("rhs_fnc",rhs,pdename_);
      ptRHSFnc_=FncReader(rhs);
      conf->get("rhs_arg_fnc",arg_rhs_,pdename_);
      SetRHSFnc=TRUE;
      directionFnc_.resize(2);
      conf->get("direction_rhs_x",directionFnc_[0],pdename_);
      conf->get("direction_rhs_y",directionFnc_[1],pdename_);
    }    
 
  }

  void Acoustic2dPDE::ComputeRHS(const Double atime)
  {
#ifdef TRACE
    (*trace) << "entering Acoustic2dPDE::ComputeRHS" << std::endl;
#endif
    Integer n;
    Integer node;
    Integer matrix_id;
    Vector<Double> coeffMass;
    Vector<Double> coeffDamp;
    std::list<Integer> nodes_hd;
    Integer i;

    coeffMass = sol_old_*a0_+sol_der1_old_*a2_+sol_der2_old_*a3_;

    algsys_->UpdateRHS(MASS,coeffMass.get());

    // damping matrix part

//     if (with_absBCs_) 
//       {

// 	coeffDamp = -sol_der1_old_-sol_der2_old_*a6_;
// 	algsys_->UpdateRHS(DAMPING,coeffDamp.get());

// 	coeffDamp = sol_old_*a1_+sol_der1_old_*a2_*a7_+sol_der2_old_*a7_*a3_;  
// 	algsys_->UpdateRHS(DAMPING,coeffDamp.get());
//       }


//     Integer level=0;          ////
//     if (SetRHSFnc){         // for standing waves

//       Double density, compress;    // get density & compressity
//       Integer i,matnum;
//       for (i=0; i<subdoms_.size(); i++)
// 	{
// 	  conf->get(subdoms_[i],matnum,"list_subdomains");
// 	  MatFile_->ReadDensityAndCompressity(density,compress,matnum,"fluid");
// 	}


//       Integer j,node1,node2;
//       Double valfnc,integrShFnc,val,multiplier;
//       Point<2> * ptCoordNodes;
//       std::vector<Double> normal;
//       normal.resize(2);
//       Vector<Integer> connecth;
//       std::vector<Elem*> edgesBC;  // vector of 1D-elements from mesh-file

//       edgesBC=ptBCs_->getEdgesBC(rhs_surfaces_[0],level);
//       valfnc = ptRHSFnc_(atime*arg_rhs_*2*3.1416); // value of fnc at the timestep

//       for (j=0; j< edgesBC.size(); j++) { // loop over surface elements
// 	ptCoordNodes=new Point<2>[2];

// 	connecth=edgesBC[j]->connect;
// 	ptgrid_->GetCoordNodesElem(connecth,ptCoordNodes,level);  
// 	calcNormal2Line(normal,ptCoordNodes[0],ptCoordNodes[1]); 

// 	multiplier=ScalarMult(normal,directionFnc_); // n * V /direction fnc/
// 	integrShFnc = edgesBC[j]->ptElem->getIntVal(ptCoordNodes);
// 	val=integrShFnc*multiplier*valfnc*density; 

// 	Double * valArr=new Double[2];
// 	valArr[0]=val;
// 	valArr[1]=val;
// 	// 2 - number of nodes in element
// 	algsys_->SetElementRHS(valArr,connecth.get(),2);

// 	delete [] valArr;
// 	delete [] ptCoordNodes;
//       }    
//     } // end of part: standing waves      

  }

}


