#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include <DataInOut/conffile.hh>
#include <DataInOut/WriteInfo.hh>
#include <multigrid.hh>

#include "assemble.hh" 

namespace CoupledField
{
  
  Assemble::Assemble(BaseSystem * algsys, Grid * aptgrid)
    :algsys_(algsys),
     ptgrid_(aptgrid),
     stiffnessMatrix_(FALSE),
     dampingMatrix_(FALSE),
     massMatrix_(FALSE),
     convectionMatrix_(FALSE),
     actlevel_(0),
     integrators_(0),
     rhsIntegrators_(0)
  {
#ifdef TRACE
    (*trace) << "entering Assemble::Assemble " << std::endl;
#endif
  }





  Assemble::~Assemble()
  {
#ifdef TRACE
    (*trace) << "entering Assemble::~Assemble " << std::endl;
#endif
    
    for (int i=0; i<integrators_.size();i++)
      for (int j=0; j<integrators_[i]->size(); j++)
	delete (*integrators_[i])[j];

    for (int i=0; i<surfintegrators_.size();i++)
      for (int j=0; j<surfintegrators_[i]->size(); j++)
	delete (*surfintegrators_[i])[j];

    for (int i=0; i<rhsIntegrators_.size();i++)
      for (int j=0; j<rhsIntegrators_[i]->size(); j++)
	delete (*rhsIntegrators_[i])[j];
  }


  Integer Assemble::SubDomIndex(const std::string & subDomName)
  {
#ifdef TRACE
    (*trace) << "entering Assemble::SubDomIndex " << std::endl;
#endif

    for (int i=0; i < subdoms_.size();i++)
      if (subDomName == subdoms_[i])
	return i;
    
  
    std::string errOut;
    errOut = "SubDomain " + subDomName + " not defined!";
    Info->Error(errOut, __FILE__, __LINE__);
  }
  
  Integer Assemble::SurfDomIndex(const std::string & surfDomName)
  {
#ifdef TRACE
    (*trace) << "entering Assemble::SurfDomIndex " << std::endl;
#endif

    for (int i=0; i < surfdoms_.size();i++)
      if (surfDomName == surfdoms_[i])
	return i;
    
  
    std::string errOut;
    errOut = "Surface-Domain " + surfDomName + " not defined!";
    Info->Error(errOut, __FILE__, __LINE__);
  }


  void Assemble::GetElemCoords(const Vector<Integer> connect, 
			       Matrix<Double> &coordMat, const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering Assemble:GetElemCoords:" << std::endl;
#endif
    
    ptgrid_->GetCoordNodesElemMat(connect, coordMat, level);
  
    //    if (deltCoords_.size() != 0)
    //      {
    //     for (Integer i=0; i<coordMat.size_row(); i++)
    //       for (Integer j=0; j<coordMat.size_col(); j++) 
    // 	coordMat(i,j) += deltCoords_[i][Mesh2PDENode_ [connect[j] - 1] - 1];
    //      }
  }

  
  // do the basic assembling stuff
  void Assemble::AssembleMatrices(const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering Assemble:AssembleMatrices" << std::endl;
#endif
    static Boolean firstTime = TRUE;
    static Boolean oneIntIsNonlin = FALSE;
    Matrix<Double> elemmat;
    const Integer nrMatrices = 5+1;  // 5 matrices, but index starts with 1!
    static Vector<Boolean> reassembleMat(nrMatrices);


    // initialize reassembling "indicator" vector
    if (firstTime)
      for (Integer actMat=0; actMat < nrMatrices; actMat++)
	reassembleMat[actMat] = FALSE;

    
     for (int actDom=0; actDom < subdoms_.size(); actDom++)
      {	
	std::vector<Elem*> elemssd;

	ptgrid_->GetElemSD(elemssd, subdoms_[actDom], level);


	for (int actEl=0; actEl< elemssd.size(); actEl++)
	  {
	    BaseFE * ptEl = elemssd[actEl]->ptElem;
	    Vector<Integer> connecth = elemssd[actEl]->connect;


	    Matrix<Double> ptCoord;
	    GetElemCoords(connecth, ptCoord, level);


	    // map connect to PDE node numbers
	    Vector<Integer> connect_PDE;
	    Mesh2PDENode(connect_PDE, connecth, *mesh2PDENode_);


	    Matrix<Double> elSol;

	    // this matrix is nonlinear and, therefore, has to be reassembled next time
	    if (oneIntIsNonlin || firstTime)
	      // fetch solution at element nodes
	      GetSolOfElement(elSol, connect_PDE);
	      
	    
	    // =========================================================================
	    //                             assemble matrices
	    // =========================================================================

	    for(Integer actInteg=0; actInteg < integrators_[actDom]->size(); actInteg++)
	      {
		IntegratorDescriptor * actDescriptor =
		  (*integrators_[actDom])[actInteg];


		// assemble only if nonlinear or first time
		if (reassembleMat[actDescriptor->DestMat()] || firstTime)
		  {
		    actDescriptor->GetIntegrator()->SetElemPtr(ptEl);
		    enum MatrixType destMat = actDescriptor->DestMat();

		    // this matrix is nonlinear and, therefore, has to be reassembled next time
		    if (actDescriptor->IsNonLin())
		      {
			oneIntIsNonlin = TRUE;
			reassembleMat[actDescriptor->DestMat()] = TRUE;
		
			actDescriptor->GetIntegrator()->SetActElemSol(elSol);
		      }

		    actDescriptor->GetIntegrator()->CalcElementMatrix(ptCoord, elemmat);
		    
		    algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), 
					      connect_PDE.size(), destMat); 
		  }		
	      }



	    // =========================================================================
	    //                             assemble RHS
	    // =========================================================================

	    for(Integer actRhsInt=0; actRhsInt < rhsIntegrators_[actDom]->size(); actRhsInt++)
	      {
		BaseIntDescriptor * actRhsID = (*rhsIntegrators_[actDom])[actRhsInt];

		actRhsID->GetIntegrator()->SetElemPtr(ptEl);

		if (actRhsID->IsNonLin())
		  actRhsID->GetIntegrator()->SetActElemSol(elSol);

		
		std::vector<Double> elemVec;
		actRhsID->GetIntegrator()->CalcElemVector(ptCoord, elemVec);
		
		// subtract internal forces on rhs from external forces 
		elemVec *= -1;


		algsys_->SetElementRHS(&elemVec[0], connect_PDE.get(), connect_PDE.size());
	      }
	    

	  }
      }

     //assemble matrices for surface integrators
     for (Integer actDom=0; actDom < surfdoms_.size(); actDom++)
      {	
	std::vector<Elem*> elemssd;
	elemssd=ptBCs_->getFacesBC(surfdoms_[actDom],level);

	for (int actEl=0; actEl< elemssd.size(); actEl++)
	  {
	    BaseFE * ptEl = elemssd[actEl]->ptElem;
	    Vector<Integer> connecth = elemssd[actEl]->connect;


	    Matrix<Double> ptCoord;
	    GetElemCoords(connecth, ptCoord, level);


	    // map connect to PDE node numbers
	    Vector<Integer> connect_PDE;
	    Mesh2PDENode(connect_PDE, connecth, *mesh2PDENode_);


	    Matrix<Double> elSol;

	    // this matrix is nonlinear and, therefore, has to be reassembled next time
	    if (oneIntIsNonlin || firstTime)
	      // fetch solution at element nodes
	      GetSolOfElement(elSol, connect_PDE);
	      
	    
	    // =========================================================================
	    //                             assemble matrices
	    // =========================================================================

	    for(Integer actInteg=0; actInteg < surfintegrators_[actDom]->size(); actInteg++)
	      {
		IntegratorDescriptor * actDescriptor =
		  (*surfintegrators_[actDom])[actInteg];


		// assemble only if nonlinear or first time
		if (reassembleMat[actDescriptor->DestMat()] || firstTime)
		  {
		    actDescriptor->GetIntegrator()->SetElemPtr(ptEl);
		    enum MatrixType destMat = actDescriptor->DestMat();

		    // this matrix is nonlinear and, therefore, has to be reassembled next time
		    if (actDescriptor->IsNonLin())
		      {
			oneIntIsNonlin = TRUE;
			reassembleMat[actDescriptor->DestMat()] = TRUE;
		
			actDescriptor->GetIntegrator()->SetActElemSol(elSol);
		      }

		    actDescriptor->GetIntegrator()->CalcElementMatrix(ptCoord, elemmat);
		    
		    algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), 
					      connect_PDE.size(), destMat); 
		  }		
	      }
	  }
      }

     firstTime = FALSE;
  }




  // do the basic assembling stuff
  void Assemble::AssembleRHS(const Integer level, const Double time)
  {
    
#ifdef TRACE
    (*trace) << "entering Assemble:AssembleRHS" << std::endl;
#endif
    AssembleRHSNodalSources(level, time);

    AssembleRHSIntegralSources(level, time);
  }
  




  // do the basic assembling stuff
  void Assemble::AssembleRHSIntegralSources(const Integer level, const Double time)
  {
#ifdef TRACE
    (*trace) << "entering Assemble:AssembleRHSIntegralSources" << std::endl;
#endif    
  
  }
  



  // do the basic assembling stuff
  void Assemble::AssembleRHSNodalSources(const Integer level, const Double time)
  {
#ifdef TRACE
    (*trace) << "entering Assemble:AssembleRHSNodalSources" << std::endl;
#endif
    
    for (int actDom=0; actDom < loadDom_.size(); actDom++)
      {
	std::string doftype = loadDom_[actDom];

	Integer dof = GetBCDof( loadDof_[actDom] );	

	std::list<Integer> nodes;
	nodes = ptBCs_->GetNodesLevel(loadDom_[actDom], level);
	
	Double val = loadVals_[actDom];

	Double val_tfunc = 1.0;
	if (ptTimeFunc_->GetmaxTimeFnc()!=0) //see, if there is any time function
	  val_tfunc = ptTimeFunc_->TimeFuncAtTime(time, level);
	
	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
	  {
	    Integer node = *p;
	    
	    val = loadVals_[actDom] * val_tfunc;
	    myCout << "val_tfunc = " << val_tfunc << myEndl;
	    myCout << "VAL = " << val << myEndl;
	    
	    algsys_->SetNodeRHS(val, (*mesh2PDENode_)[node-1], dof);	
	  }
      }
  }
  



Integer Assemble::
GetBCDof(const std::string dofString)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::GetBCDof " << std::endl;
#endif

  if (dofString == "ux")
    return 1;
  else
    if (dofString == "uy")
      return 2;
    else
      if (dofString == "uz")
	return 3;
      else
	Error("The direction mentioned in the config-file is not implemented! ",__FILE__,__LINE__);
}




  void Assemble::Mesh2PDENode(Vector<Integer> & PDENodes, 
			      const Vector<Integer> & MeshNodes,
			      const std::vector<Integer> & Mesh2PDENode)
  {
#ifdef TRACE
    (*trace) << "entering Assemble::Mesh2PDENode " << std::endl;
#endif
    
    PDENodes.Resize(MeshNodes.size());
    
    for (Integer i=0; i<MeshNodes.size(); i++) 
      PDENodes[i] = Mesh2PDENode[MeshNodes[i]-1];
  }  





  void Assemble::InitMatrices()
  {
    // Initialize matrices in order to get BCs correct
    algsys_->InitMatrix(SYSTEM);

    if (stiffnessMatrix_)
      algsys_->InitMatrix(STIFFNESS);
    
    if (dampingMatrix_)
      algsys_->InitMatrix(DAMPING);

    if (convectionMatrix_)
      algsys_->InitMatrix(CONVECTION);
    
    if (massMatrix_)
      algsys_->InitMatrix(MASS); 
  }
 
 

  void Assemble::CreateMatrices()
  {
    const Integer numconstraints = 0;  // currently not handled
    
    Integer matrixsystype[5];
    
    matrixsystype[0] = SYSTEM;      // memory for the system matrix
    if (stiffnessMatrix_  == 1) matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
    if (dampingMatrix_    == 1) matrixsystype[2] = DAMPING;     // memory for the damping matrix
    if (convectionMatrix_ == 1) matrixsystype[3] = CONVECTION;  // memory for the convection matrix
    if (massMatrix_       == 1) matrixsystype[4] = MASS;        // memory for the mass matrix
    
    //put to algebraic system
    algsys_->CreateMatrix(matrixsystype, matrixType_, graphType_, 
			  dofsPerNode_,numDirichletBCs_, numconstraints);
  }
  

  void Assemble::SetGeneralParams(const std::string & pdename, 
				  const Integer dofsPerNode,
				  const Integer numPDENodes, 
				  const std::vector<std::string> subdoms,
				  const std::vector<std::string> surfdoms)
  {
    pdename_     = pdename;
    dofsPerNode_ = dofsPerNode;
    numPDENodes_ = numPDENodes;
    subdoms_     = subdoms;
    surfdoms_    = surfdoms;



    // read load values =========================================
    conf->ifgetliststr("loads", loadDom_, pdename_);
    conf->ifgetliststr("loadDof", loadDof_, pdename_);
    

    if (loadDom_.size() != loadDof_.size())   //check for load data
      Error("Inconsistent definition of loads");
    
    loadVals_.resize(loadDom_.size());
    
    for(int i=0; i < loadDom_.size(); i++)
      conf->get(loadDom_[i], loadVals_[i], pdename_, "bc_conditions","loads");
    



    // for every domain, we need an own integrator list ==========
    integrators_.resize(subdoms_.size());
    surfintegrators_.resize(surfdoms_.size());
    rhsIntegrators_.resize(subdoms_.size());
    for (int i=0; i<subdoms_.size();i++)
      {
	integrators_[i] = new std::vector<IntegratorDescriptor *>;
	surfintegrators_[i] = new std::vector<IntegratorDescriptor *>;
	rhsIntegrators_[i] = new std::vector<BaseIntDescriptor *>;
      }
  }
  

  //  void Assemble::SetupMatrixGraph(Integer numeq, Integer graphType)
  void Assemble::SetupMatrixGraph(Integer numeq)
  {
#ifdef TRACE
    (*trace) << "entering Assemble::SetupMatrixGraph" << std::endl;
#endif
    
  //initialize matrix graph
  Integer estimated_maxneighbors=20;
  algsys_->CreateGraph(numeq, graphType_,estimated_maxneighbors);

  // set the graph - connectivity matrix
  BaseFE * ptElem; 
  Integer nsub, iel, fe_type;
  Vector<Integer> connecth;

  for (nsub=0; nsub<subdoms_.size(); nsub++)
    {
      std::vector<Elem*> elemssd;
      ptgrid_->GetElemSD(elemssd,subdoms_[nsub],actlevel_);

      for (iel=0; iel < elemssd.size(); iel++)
	{  
	  ptElem=elemssd[iel]->ptElem;
	  
	  //Map Mesh Node numbers to PDE node numbers
	  Mesh2PDENode(connecth,elemssd[iel]->connect, *mesh2PDENode_);

	  fe_type=elemssd[iel]->ptElem->feType();
	  algsys_->SetElementPos(connecth.get(),connecth.size(),fe_type);
	}
    }

}



  Integer Assemble::GetNrBCDof(const std::string & dofStartString)
  {
#ifdef TRACE
    (*trace) << "entering Analysis::GetNrBCDof" << std::endl;
#endif
    
    Integer nrActDof;
    
    if (dofStartString == "ux")
      nrActDof = 1;
    else 
      if (dofStartString == "uy")
	nrActDof = 2;
      else
	if (dofsPerNode_ == 3)
	  if (dofStartString == "uz")
	    nrActDof = 3;
	  else
	    Error("Unknown dof-type in homog. BC; substring must start with ux, uy or uz!!",__FILE__,__LINE__);
	else
	  Error("Unknown dof-type in homog. BC; substring must start with ux or uy!!",__FILE__,__LINE__);
    
    return nrActDof;
  }
  


void Assemble::
GetSolOfElement( Matrix<Double>& elDisp, Vector<Integer>& connect_PDE)
{
#ifdef TRACE
    (*trace) << "entering Assemble::GetSolOfElement" << std::endl;
#endif

  // displacement of element nodes
  elDisp.Resize(dofsPerNode_, connect_PDE.size());

  for (Integer dim=0; dim<dofsPerNode_; dim++)
    for(Integer actNode=0; actNode<connect_PDE.size(); actNode++)
      elDisp[dim][actNode] = (*sol_)[dim][connect_PDE[actNode]-1];
}





/// define integrators
void Assemble::AddRhsIntegrator(BaseForm * integrator, const std::string & subDomName, 
				const Integer nonLin)
{
#ifdef TRACE
  (*trace) << "entering Assemble::AddRhsIntegrator " << std::endl;
#endif

    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsIntegrators_[SubDomIndex(subDomName)]->push_back(actRhsID);
  }


  // ==========================================================
  // STATIC ANALYSIS
  // ==========================================================


  StaticAssemble::StaticAssemble(BaseSystem * algsys, Grid * agrid)
    :Assemble(algsys, agrid)
  {
    graphType_    = NODEGRAPH; 
  }
  


  /// define integrators
  void StaticAssemble::AddIntegrator(BaseForm * integrator, const std::string & subDomName,
					const enum MatrixType destinationMatrix, const Integer nonLin)
  {
#ifdef TRACE
    (*trace) << "entering StaticAssemble::AddIntegrator " << std::endl;
#endif
    MatrixType actMatType = destinationMatrix;
    
    if (actMatType == STIFFNESS)
      actMatType = SYSTEM;

    if (actMatType !=  SYSTEM)
      return;

    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator, actMatType, nonLin);
    integrators_[SubDomIndex(subDomName)]->push_back(actID);
  }

  /// define integrators
  void StaticAssemble::AddSurfIntegrator(BaseForm * integrator, const std::string & subDomName,
					const enum MatrixType destinationMatrix, const Integer nonLin)
  {
#ifdef TRACE
    (*trace) << "entering StaticAssemble::AddSurfIntegrator " << std::endl;
#endif
    MatrixType actMatType = destinationMatrix;
    
    if (actMatType == STIFFNESS)
      actMatType = SYSTEM;

    if (actMatType !=  SYSTEM)
      return;

    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator, actMatType, nonLin);
    surfintegrators_[SurfDomIndex(subDomName)]->push_back(actID);
  }




  /// define integrators
  void StaticAssemble::AddIntegrator(IntegratorDescriptor * actID, const std::string & subDomName)
  {
#ifdef TRACE
    (*trace) << "entering StaticAssemble::AddIntegrator " << std::endl;
#endif
    
    
    if (actID->DestMat() == STIFFNESS)
      actID->SetDestMat(SYSTEM);

    if (actID->DestMat() !=  SYSTEM)
      return;

    integrators_[SubDomIndex(subDomName)]->push_back(actID);
  }
    


  // ==========================================================
  // TRANSIENT ANALYSIS
  // ==========================================================


  TransientAssemble::TransientAssemble(BaseSystem * algsys, Grid * agrid)
    :Assemble(algsys, agrid)
  {
     graphType_    = NODEGRAPH; 

    stiffnessMatrix_  = TRUE;
    massMatrix_       = TRUE;
  }



  /// define integrators
  void TransientAssemble::AddIntegrator(BaseForm * integrator, const std::string & subDomName,
					const enum MatrixType destinationMatrix, const Integer nonLin)
  {
#ifdef TRACE
    (*trace) << "entering TransientAssemble::AddIntegrator " << std::endl;
#endif
    if (destinationMatrix == SYSTEM)
      Info->Error("In transient assembling, no SYSTEM matrix may be defined directly", __FILE__, __LINE__);

    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator, destinationMatrix, nonLin);
    integrators_[SubDomIndex(subDomName)]->push_back(actID);
  }

  /// define integrators
  void TransientAssemble::AddSurfIntegrator(BaseForm * integrator, const std::string & subDomName,
					const enum MatrixType destinationMatrix, const Integer nonLin)
  {
#ifdef TRACE
    (*trace) << "entering TransientAssemble::AddSurfIntegrator " << std::endl;
#endif
    if (destinationMatrix == SYSTEM)
      Info->Error("In transient assembling, no SYSTEM matrix may be defined directly", __FILE__, __LINE__);

    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator, destinationMatrix, nonLin);
    surfintegrators_[SurfDomIndex(subDomName)]->push_back(actID);
  }

  /// define integrators
  void TransientAssemble::AddIntegrator(IntegratorDescriptor * actID, const std::string & subDomName)
  {
#ifdef TRACE
    (*trace) << "entering TransientAssemble::AddIntegrator " << std::endl;
#endif
    if (actID->DestMat() == SYSTEM)
      Info->Error("In transient assembling, no SYSTEM matrix may be defined directly", __FILE__, __LINE__);

    integrators_[SubDomIndex(subDomName)]->push_back(actID);
  }




  // ==========================================================
  // RHS INTEGRATOR DESCRIPTOR
  // ==========================================================

BaseIntDescriptor::BaseIntDescriptor()
  :integrator(NULL),
   nonLin(FALSE)
{
#ifdef TRACE
    (*trace) << "entering BaseIntDescriptor::BaseIntDescriptor" << std::endl;
#endif
  }


  BaseIntDescriptor::~BaseIntDescriptor()
  {
#ifdef TRACE
    (*trace) << "entering BaseIntDescriptor::~BaseIntDescriptor" << std::endl;
#endif

    if (integrator)
      delete integrator;
  }
  

  /// define integrators
  BaseIntDescriptor::BaseIntDescriptor(BaseForm * aIntegrator, const Boolean aNonLin)
    :integrator(aIntegrator),
     nonLin(aNonLin)
  {
#ifdef TRACE
    (*trace) << "entering BaseIntDescriptor::BaseIntDescriptor " << std::endl;
#endif

  }
  



  // ==========================================================
  // INTEGRATOR DESCRIPTOR
  // ==========================================================


  IntegratorDescriptor::IntegratorDescriptor()
    :BaseIntDescriptor(),
     destinationMatrix(SYSTEM),
     secondaryMatrix(NOTYPE),
     secMatFac(0.0)
{
}
  

  /// define integrators
  IntegratorDescriptor::IntegratorDescriptor(BaseForm * aIntegrator, 
					     const enum MatrixType aDestMat, const Boolean aNonLin)
    :BaseIntDescriptor(aIntegrator, aNonLin),
     destinationMatrix(aDestMat)

  {
#ifdef TRACE
    (*trace) << "entering IntegratorDescriptor::IntegratorDescriptor " << std::endl;
#endif
  }
  





  IntegratorDescriptor::~IntegratorDescriptor()
  {
    if (integrator)
      delete integrator;
  }
  


} // end namespace CoupledField
