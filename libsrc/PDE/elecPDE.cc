#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>

#include "elecPDE.hh"
#include <DataInOut/Unverg/outUnverg.hh>
#include <Forms/forms_header.hh>
#include <Forms/elecfieldop.hh>
#include <Forms/elecforceop.hh>
#include <Estimator/spaceerror.hh>
 
namespace CoupledField
{

ElecPDE::ElecPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
:BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::ElecPDE " << std::endl;
#endif

  dofspernode_=1;
  E_ = new Vector<Double>[ptgrid_->GetDim()];
  Force_ = new Vector<Double>[ptgrid_->GetDim()];

  SetMatrixFactors();

  pdename_    = "electrostatic";
  pdematerialclass_ = "piezo";
  
  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);

  AssignPDENodeNumbers();
  size_ = NumPDENodes_;

  NumElems_ = ptgrid_->GetMaxnumElem(actlevel_, subdoms_); 
  sol_.Resize(size_);
  sol_.Init(0);

}


void ElecPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 0.0;  // factor for stiffness matrix
  matrix_factor_[1] = 0.0;  // factor for damping matrix
  matrix_factor_[2] = 0.0;  // factor for convection matrix
  matrix_factor_[3] = 0.0;  // factor for mass matrix
}


void ElecPDE::DiscreteParamsPDE()
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::DiscreteParamsPDE" << std::endl;
#endif

  MatrixType_   = RSCALAR;
  GraphType_    = NODEGRAPH; 
  SystemMatrix_ = TRUE;

}

void ElecPDE::SolveStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::SolveStepStatic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 1;

  Double * ptsol;

  //compute and assemble element matrices
  SetupMatrices(level);

  //account for bcs
  SetBCs(level,update,0);

  algsys_->CalcPrecond(job);
  algsys_->Solve();

  ptsol = algsys_->GetSolutionVal();

  // save solution
  Vector<Double> transsol(NumPDENodes_, ptsol);
  sol_=transsol;

#ifdef DEBUG
  //  std::string matFileName = "solMat";
  //  OutFile_->WriteSolMatrix(ptgrid_, level, sol_, matFileName);
#endif
}


void ElecPDE::PostProcess(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::PostProcess" << std::endl;
#endif  

  ShortInt Dim = ptgrid_->GetDim();
  ElecFieldOp * FieldOp = new ElecFieldOp(ptgrid_, &Mesh2PDENode_, &sol_, level);
  ElecForceOp * ForceOp = new ElecForceOp(ptgrid_, &Mesh2PDENode_, &sol_, level);
 
  // ######### Calculation of the electric field #########

  std::vector<Double> LCoord;
  LCoord.resize(2);
  LCoord[0] = 0;
  LCoord[1] = 0;

  std::vector<Elem*> elemssd;
  Integer counterElems=0;
  Vector<Double> TempE, TempFP;
  
  // Resize solution arrays
  for( Integer i=0; i<Dim; i++) 
    {
      E_[i].Resize(NumElems_);
      Force_[i].Resize(NumElems_);
    }
  
  // loop over all subdomains
  for (Integer isd=0; isd<subdoms_.size(); isd++)
    {
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid_->GetElemSD(elemssd,subdoms_[isd],level);
      
      // loop over elements of subdomain
      for (Integer iel=0; iel< elemssd.size(); iel++,counterElems++)
	{
 	  FieldOp->CalcElemElecField( TempE, elemssd[iel], LCoord);
	  for (Integer k=0; k<Dim; k++)
	    
	    E_[k][counterElems] = TempE[k];
	}
    }

//   // ########### Calculation of Electric Force ########

//   std::vector<std::string> force_surfaces;
//   conf->getliststr("force_surfaces",force_surfaces,pdename_);

//  //For the subdomain next to boundary of obstacle (Next2Surf)
//   //std::vector<Elem*> Next2Surf;  // vector of 2D-elements from mesh-file      
//   //BaseFE * ptElNext2Surf;
//   std::vector<Elem*> belongSE, InterfaceSurf, Next2Surf;
//   std::vector<ShortInt> IsBoundaryNode;
//   std::vector<Integer> map;
  
//   InterfaceSurf=ptBCs_->getEdgesBC(force_surfaces[0],level);
  
//   // TESTING Routine for Interface Surface
//   // for( Integer i=0; i<InterfaceSurf.size(); i++)
//   //  std::cerr << "InterfaceSurf[" << i<< "].connect = " << InterfaceSurf[i]->connect << std::endl;
  
//   //std::cerr << "force_surfaces[1] = " << force_surfaces[1] << std::endl;


//   Next2Surf=ptBCs_->getNeighElemsForSurfaces(force_surfaces[1],level);

//   //for( Integer i=0; i<Next2Surf.size(); i++)
//   //  std::cerr << "Next2Surf[" << i<< "].connect = " << Next2Surf[i]->connect << std::endl;


//   //ptgrid_->DefineBelonging4Elems(InterfaceSurf,Next2Surf,belongSE);
//   ptgrid_->GetInterfaceNeighbours(InterfaceSurf, Next2Surf, belongSE);

//   //std::cerr << "belongSE.size() = " << belongSE.size() << std::endl;
//   //for( Integer i=0; i<belongSE.size(); i++)
//   //  std::cerr << "belongSE[" << i<< "].connect = " << belongSE[i]->connect << std::endl;

  
//   // Get Material Parameter
//   Double epsilon;

//   for (Integer i=0; i<subdoms_.size(); i++)
//     if (force_surfaces[1] == subdoms_[i]) 
//       epsilon = materialData_[i].GetPermittivity(2,2); 

//   // Get Node numbers for boundary interface
//   ptgrid_->CalcNumberOfNodesInPatch(InterfaceSurf, map);
  

//   // loop over all elements of moving surface boundary
//   //std::cerr << "belongeSE.size() = " << belongSE.size() << std::endl;
//   for (Integer isd=0; isd<belongSE.size(); isd++)
//     {

//       IsBoundaryNode.resize(belongSE[isd]->connect.size(), 0);
      
      
//       // Determine Boundary Nodes
//       for (Integer inode=0; inode<IsBoundaryNode.size(); inode++)
// 	{
// 	  for (Integer imap=0; imap <map.size(); imap++)
// 	    {
// 	      if (belongSE[isd]->connect[inode] == map[imap] )
// 		{
// 		  IsBoundaryNode[inode] = 1;
// 		  break;
// 		}
// 	    }
// 	}
      
//       // std::cerr << "belongeSE[" << isd << "].connect = " << belongSE[isd]->connect << std::endl;

//       //for (Integer d=0; d<IsBoundaryNode.size(); d++)
//       //std::cerr << IsBoundaryNode[d] << " ";
//       //std::cerr <<  std::endl;

//       ForceOp->CalcElemElecForce( TempFP, belongSE[isd], epsilon, IsBoundaryNode);
      
//       for (Integer k=0; k<Dim; k++)
// 	Force_[k][isd] = TempFP[k];
      
//     }

}


void ElecPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::SetupMatrices" << std::endl;
#endif
  
  Matrix<Double> elemmat;  
  Matrix<Double> ptCoord;
  BaseFE         * ptElem;

  if (InfoPrint)
    (*infofile) << " ------------------------- Element matrices --------------- " << std::endl;

  Vector<Integer> connecth, connect_PDE;  
  Integer i, j;

  for (i=0; i<subdoms_.size(); i++)
    {
      //reads eps33 (matrix notation starts with 0)
      Double eps33 = materialData_[i].GetPermittivity(2,2);

      std::vector<Elem*> elemssd;
   
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j < elemssd.size(); j++)
	{  
	  ptElem=elemssd[j]->ptElem;

	  BaseForm * bilinear_stiff = new LaplaceInt(ptElem, eps33);

	  connecth=elemssd[j]->connect;
	  
	  ptgrid_->GetCoordNodesElemMat(connecth, ptCoord, level);

	  // CHANGE connecth
	  Mesh2PDENode(connect_PDE,connecth);

	  // stiffness part
	  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);
	  
#ifdef DEBUG
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif

	  algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connecth.size(), SYSTEM);
	  
	  delete bilinear_stiff;	  

	}  
      
    }

}


void ElecPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::WriteResultsInFile" << std::endl;
#endif

  Integer step=0;
  Double time=0;

  ////////////////////////
  //// CHANGE   //////////
  ////////////////////////

  Vector<Double> Potentialh, *Eh, *Fh;
  
  TransformNodeSolution(Potentialh,sol_);
  Eh = new Vector<Double>[ptgrid_->GetDim()];
  Fh = new Vector<Double>[ptgrid_->GetDim()];
  for (Integer i=0; i<ptgrid_->GetDim(); i++) 
    {
      TransformElemSolution(Eh[i],E_[i]);
      TransformElemSolution(Fh[i],Force_[i]);
    }

  // Transform solution from PDE to Mesh nodes

  // write results
  if (OutFile_->IsGMV())
    {
    OutFile_->WriteSolution(Potentialh,step,time,"electric_potential");

    // Write Out Vector Data
    for (ShortInt i=0; i<ptgrid_->GetDim(); i++) 
      {
	std::ostringstream el_fieldname, f_fieldname;
	el_fieldname << "Efield" << i;
	//	OutFile_->WriteDataOnCell(Eh[i],step,time, el_fieldname.str());
	f_fieldname << "Eforce" << i;
	//	OutFile_->WriteDataOnCell(Fh[i],step,time, f_fieldname.str());
      }
    }
  else
    OutFile_->WriteSolution(Potentialh,step,time,"electric_potential");

}

ElecPDE::~ElecPDE()
{
  if (E_) delete[] E_;
  if (Force_) delete[] Force_;
}

} // end of namespace


