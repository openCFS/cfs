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
 
#include <AlgebraicSystem/LinAlg/linsystem.hh>

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
  
  SetMatrixFactors();

  pdename_    = "electrostatic";
  pdematerialclass_ = "piezo";
  
  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);

  Reset();
  
  SolverCFS_ = FALSE;
}

void ElecPDE::Reset()
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::SetMatrixFactors" << std::endl;
#endif
    
  AssignPDENodeNumbers(Mesh2PDENode_,PDE2MeshNode_,subdoms_);
  NumPDENodes_=PDE2MeshNode_.size();

  deltCoords_.reshape(Dim_,NumPDENodes_);
  sol_.reshape(dofspernode_,NumPDENodes_);
  sol_.init();

  NumElems_=ptgrid_->GetMaxnumElem(actlevel_,subdoms_);

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


void ElecPDE::InitCoupling(PDECoupling * Coupling)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::InitCoupling" << std::endl;
#endif

  ptCoupling_ = Coupling;


  // Initialization of coupling helper arrays
  std::string quantity;
  std::vector<Integer> * couplingnodes;
  std::vector<Elem*> interface_tmp;
  std::vector<std::vector<ShortInt> > isBoundaryNode_tmp;
  //std::vector<Integer> numBoundaryNodes_tmp;
  std::vector<std::vector<Integer> > elemNodeToCouplingNode_tmp;
  Array<Double> * values;

  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      // Initialize arrays for coupling surface elements
      if (ptCoupling_->GetOutputQuantity(i) == "elecforce")
	{
	  ptCoupling_->GetOutputNodes(i, couplingnodes);
	  if (couplingnodes == 0)
	    std::cerr << "Couplingnodes = 0!!!!" << std::endl;
	  
	  ptgrid_->GetInterfaceNeighbours(*couplingnodes, subdoms_, interface_tmp, actlevel_);
	  
	  F_Interface_.push_back(interface_tmp);

	  // Intialize the memory
	  ptCoupling_->SetOutputDim(i, Dim_);
	 
	  // Assign arrays for element boundary nodes
	  isBoundaryNode_tmp.clear();
	  isBoundaryNode_tmp.resize(interface_tmp.size());
	  elemNodeToCouplingNode_tmp.clear();
	  elemNodeToCouplingNode_tmp.resize(interface_tmp.size());
	  //numBoundaryNodes_tmp.clear();
	  //numBoundaryNodes_tmp.resize(interface_tmp.size());
	 
	  
	  for (Integer ielem=0; ielem<interface_tmp.size(); ielem++)
	    {
	      //std::cerr << "Interface[" << ielem << "]: " << interface_tmp[ielem]->ElemNum << std::endl;

	      isBoundaryNode_tmp[ielem].clear();
	      isBoundaryNode_tmp[ielem].resize(interface_tmp[ielem]->connect.size());
	      elemNodeToCouplingNode_tmp[ielem].clear();
	      elemNodeToCouplingNode_tmp[ielem].resize(interface_tmp[ielem]->connect.size());

	      // Determine Boundary Nodes
	      for (Integer ielemnode=0; ielemnode<isBoundaryNode_tmp[ielem].size(); ielemnode++)
		{
		  for (Integer inodes=0; inodes<(*couplingnodes).size(); inodes++)
		    {
		      if (interface_tmp[ielem]->connect[ielemnode] == (*couplingnodes)[inodes] )
			{
			  isBoundaryNode_tmp[ielem][ielemnode] = 1;
			  elemNodeToCouplingNode_tmp[ielem][ielemnode] = inodes;
			  //numBoundaryNodes_tmp[ielem] +=1;
			  break;
			} // end if
		    } // end for (inodes)
		} // end for (ielemnodes)
	    } // end for (ielems)

	  isBoundaryNode_.push_back(isBoundaryNode_tmp);
	  elemNodeToCouplingNode_.push_back(elemNodeToCouplingNode_tmp);
	  //numBoundaryNodes_.push_back(numBoundaryNodes_tmp);

	    
// 	  for (Integer k=0; k<isBoundaryNode_tmp.size(); k++)
// 	    {
// 	      //std::cerr << "Element " << interface_tmp[k]->ElemNum << " nodes: " << interface_tmp[k]->connect << std::endl;
// 	      for (Integer l=0; l<4; l++)
// 		std::cerr << isBoundaryNode_tmp[k][l] << " ";

// 	      std::cerr << std::endl << std::endl;
// 	    }
	  
	} // end if
            
    } // end for (i)
}
  

void ElecPDE::SolveStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::SolveStepStatic" << std::endl;
#endif

  Integer job = 1;

  Double * ptsol;

  //compute and assemble element matrices
  SetupMatrices(level);

  //account for bcs
  SetBCs(level,updateBCs_,0);

  updateBCs_ = 0;
  
  algsys_->CalcPrecond(job);
  algsys_->Solve();

  ptsol = algsys_->GetSolutionVal();

  // save solution
  Integer k=0;
  
  for (Integer i=0; i<NumPDENodes_; i++)
    for (Integer dim=0; dim<dofspernode_; dim++)
      sol_[dim][i] = ptsol[k++];

  // Initalize RHS and Solution vector
  algsys_->InitRHS();
  algsys_->InitSol();  

  // Initalize matrices, if PDE couples via COORD or MAT

  if (InitMatrices_ == TRUE)
    {
      // after merging with FRED's Code, unqote following line
      // and delete the rest of the lines
      // InitMatrices();

      // Temporarily
      Integer matrixsystype[5];    
      if (SystemMatrix_     == 1) matrixsystype[0] = SYSTEM;      // memory for the system matrix
      if (StiffnessMatrix_  == 1) matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
      if (DampingMatrix_    == 1) matrixsystype[2] = DAMPING;     // memory for the damping matrix
      if (ConvectionMatrix_ == 1) matrixsystype[3] = CONVECTION;  // memory for the convection matrix
      if (MassMatrix_       == 1) matrixsystype[4] = MASS;        // memory for the mass matrix
          
      for (Integer i=0;i<5;i++)
	{
	  if (matrixsystype[i] !=0)
	    algsys_->InitMatrix(i+1);
	}
    }

  if (flags->CalcErrorMap_)
    {
      Double         totalErr;
      Array<Double>  Sol_Mesh;
      Vector<Double> solVec;
      
      ptError_=new SpaceErrorEstimator();

      ptError_->Init(this);
      
      TransformNodeSolution(Sol_Mesh,sol_,PDE2MeshNode_);

      solVec.Resize(sol_.size());
      int i;
      for (i=0; i<sol_.size(); i++)
	solVec[i]=sol_[0][i];
      //      sol_.toVector(solVec,1);

      ptError_->CalcErrorMap(solVec,subdoms_,ptgrid_,errorMap_,totalErr,level);
      
      *infofile << " total error of calculation:: " << totalErr << std::endl;
      *data << errorMap_ << std::endl;
    }

}


void ElecPDE::CalcOutputCoupling()
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::CalcOutputCoupling" << std::endl;
#endif

  std::string quantity;
  std::vector<Integer> * couplingnodes;
  Array<Double> * values;
  Integer forcesCount = 0;
  
  // loop over all output coupling quantities
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      quantity = ptCoupling_->GetOutputQuantity(i);
      
      switch(ptCoupling_->GetOutputType(i))
	{
	  
	case NODE:
	  
	  ptCoupling_->GetOutputNodes(i, couplingnodes);
	  ptCoupling_->GetOutputValues(i, values);
	  
	  if (quantity == "elecpotential")
	    NodeSolutionToCoupling(*values, *couplingnodes);
	    
	  if (quantity == "elecforce")
	    {
	      CalcNodeForce(*values, 
			    *couplingnodes, 
			    F_Interface_[forcesCount], 
			    isBoundaryNode_[forcesCount], 
			    elemNodeToCouplingNode_ [forcesCount]);
	      
	      forcesCount++;
	    }
	  
	  break;

	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}

    }

}


void ElecPDE::PostProcess(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::PostProcess" << std::endl;
#endif  

   ElecFieldOp * FieldOp = new ElecFieldOp(ptgrid_, this, &Mesh2PDENode_, &sol_[0], level);

  // ------ Calculation of the electric field ------

  std::vector<Double> LCoord;
  LCoord.resize(Dim_);
  LCoord[0] = 0;
  LCoord[1] = 0;

  std::vector<Elem*> elemssd;
  Integer counterElems=0;
  Vector<Double> TempE;
  
  // Resize solution arrays
  E_.reshape(Dim_,NumElems_);
  E_.init();
  
  
  // loop over all subdomains
  for (Integer isd=0; isd<subdoms_.size(); isd++)
    {
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid_->GetElemSD(elemssd,subdoms_[isd],level);
      
      // loop over elements of subdomain
      for (Integer iel=0; iel< elemssd.size(); iel++,counterElems++)
	{
 	  FieldOp->CalcElemElecField( TempE, elemssd[iel], LCoord);
	  
	  E_.setValuesRow(TempE,counterElems);
	}
    }
   delete FieldOp;
}


void ElecPDE::CalcNodeForce(Array<Double> & force, 
			    std::vector<Integer> & nodes, 
			    std::vector<Elem*> & elems,
			    std::vector<std::vector<ShortInt> > & isBoundaryNode,
			    std::vector<std::vector<Integer> > & elemNodeToCouplingNode)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::CalcNodeForce" << std::endl;
#endif  
  
  
  ElecForceOp * ForceOp = new ElecForceOp(ptgrid_, this, &Mesh2PDENode_, &sol_[0], actlevel_);
   
  Array<Double> force_temp;
  force.reshape(Dim_, nodes.size());
  force.init();
  
  for (Integer ielem=0; ielem<elems.size(); ielem++)
    {
      // Get Material Parameter
      Double epsilon;
      
      for (Integer i=0; i<subdoms_.size(); i++)
	{
	  if (elems[i]->namesd == subdoms_[i]) 
	  epsilon = materialData_[i].GetPermittivity(2,2); 
	}
      
      // Only for testing
      //epsilon = 8.854e-12;

      ForceOp->CalcElemElecForce( force_temp, elems[ielem], epsilon, isBoundaryNode[ielem]);
      

      // Add the element force to the according coupling node
      for (Integer ielemnode=0; ielemnode<elems[ielem]->connect.size(); ielemnode++)
	for( Integer idim=0; idim<Dim_; idim++)
	  force[idim][elemNodeToCouplingNode[ielem][ielemnode]] += force_temp[idim][ielemnode];

    }

  Vector<Double> sum;
  sum.Resize(Dim_);
  
  for (Integer i=0; i<nodes.size(); i++) 
    {
      //std::cerr << "Node[" << nodes[i] << "] = ";
    for (Integer dim=0; dim<Dim_; dim++)
      {
	sum[dim] += force[dim][i];
	//std::cerr << force[dim][i] << " , ";
      }
    //std::cerr << std::endl;
    }
  delete ForceOp;


  //std::cerr << "Sum of E-Force:" << std::endl << sum << std::endl;
  
  if (! InfoPrint)
    return;

    // write information in .info-file
 
  (*infofile) << "Sum of electrostatic force:" << std::endl;
  (*infofile) << sum << std::endl;

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
	  
	  GetElemCoords(connecth, ptCoord, level);

	  // Mape Mesh to PDE node numbers
	  Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_);

	  // stiffness part
	  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);
	  
#ifdef DEBUG
	  (*debug) << "Stiffnessmatrix, ElementNumber  " << j << std::endl;
	  (*debug) << elemmat << std::endl;
#endif

	  if (SolverCFS_)
	    {
	      Integer elsize=connecth.size();
	      Integer irow,icln;
	      Integer ii, iii;
	      for (ii=0; ii<elsize; ii++)
		for (iii=0; iii<elsize; iii++)
		  {
		    irow=connecth[ii]-1;
		    icln=connecth[iii]-1;
		    sysmat_.Add(irow,icln,elemmat[ii][iii]);
		  } 
	    }
	  
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
  ShortInt Dim = ptgrid_->GetDim();

  Array<Double> E_Mesh, Force_Mesh, Sol_Mesh;
  
  // transform solution vector for electric potential
  TransformNodeSolution(Sol_Mesh,sol_,PDE2MeshNode_);

  TransformElemSolution(E_Mesh,E_,subdoms_);

  // CHANGE F_Interface_
  // TransformElemSolution(Force_Mesh,Force_,F_Interface_[0]);
   
  // write results
  if (OutFile_->IsGMV())
    {
      // write electric potential
      OutFile_->WriteNodeSolution(Sol_Mesh,step,time,"E-Potential");
      
      // Write Out Vector Data
      OutFile_->WriteElemSolution(E_Mesh,step,time,"E-Field");
      //OutFile_->WriteElemSolution(Force_Mesh,step,time,"E-Force");
    }
  else
    {
      // write electric potential
      OutFile_->WriteNodeSolution(Sol_Mesh,step,time,"electric potential");
      
      // Write Out Vector Data
       OutFile_->WriteElemSolution(E_Mesh,step,time,"electric field");
      //OutFile_->WriteElemSolution(Force_Mesh,step,time,"mag. volume force");
    }

    if (flags->CalcErrorMap_)
    {   
      OutFile_->WriteElemSolution(errorMap_,step,time,"relERR-E-Potential"); 
    }
}

Boolean ElecPDE::HasOutput(std::string output)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::HasOutput" << std::endl;
#endif
  
  if (output == "elecforce")
    return TRUE;

  if (output == "elecpotential")
    return TRUE;

  if (output == "elecfield")
    return TRUE;

  return FALSE;

}



ElecPDE::~ElecPDE()
{
  ;
}

} // end of namespace


