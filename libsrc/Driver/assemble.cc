#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include <Domain/elem.hh>
#include <DataInOut/ParamHandling/ConfFile.hh>
#include <DataInOut/ParamHandling/BaseParamHandler.hh>
#include <DataInOut/WriteInfo.hh>


#ifdef USE_OLAS
#include <olas.hh>
#else
#include <multigrid.hh>
#endif

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
     rhsIntegrators_(0),
     rhsSrcIntegrators_(0)
  {
    ENTER_FCN( "Assemble::Assemble" );
    firstTime_ = TRUE;
    oneIntIsNonlin_ = FALSE;
    nrMatrices_ = 5+1;  // 5 matrices, but index starts with 1!
    reassembleMat_.Resize(nrMatrices_);
    nonLinGeo = FALSE;
    deltaCoords_ = NULL;

#ifdef USE_OLAS
    olasParams_ = algsys_->GetOLASParams();
    olasReport_ = algsys_->GetOLASReport(); 
#endif

  }





  Assemble::~Assemble()
  {
    ENTER_FCN( "Assemble::~Assemble" );
    
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
    ENTER_FCN( "Assemble::SubDomIndex" );

    for (int i=0; i < subdoms_.size();i++)
      {
	if (subDomName == subdoms_[i]) return i;
      }
  
    std::string errOut;
    errOut = "SubDomain " + subDomName + " not defined!";
    Info->Error(errOut, __FILE__, __LINE__);

    return 0;
  }
  
  Integer Assemble::SurfDomIndex(const std::string & surfDomName)
  {
    ENTER_FCN( "Assemble::SurfDomIndex" );

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
    ENTER_FCN( "Assemble:GetElemCoords" );
    
    ptgrid_->GetCoordNodesElemMat(connect, coordMat, level);
  
    if (nonLinGeo)
      {
	if (deltaCoords_ == NULL)
	  Error("ElecPDE: set input_coupling_terms = smoothdisplacement or nonlin = no");

	Vector<Integer> connect_PDE;
	Mesh2PDENode(connect_PDE, connect, *mesh2PDENode_);
	Double val;
	for (Integer i=0; i<coordMat.GetSizeRow(); i++)
	  for (Integer j=0; j<coordMat.GetSizeCol(); j++) 
	    {
	      val = (*deltaCoords_)[i][connect_PDE[j] - 1];
	      coordMat(i,j) += val;
	    }
      }
  }

  
  // do the basic assembling stuff
  void Assemble::AssembleMatrices(const Integer level)
  {

    ENTER_FCN( "Assemble:AssembleMatrices" );
    
    Matrix<Double> elemmat;

    // initialize reassembling "indicator" vector
    if (firstTime_)
      for (Integer actMat=0; actMat < nrMatrices_; actMat++)
	reassembleMat_[actMat] = FALSE;

    
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
	    if (oneIntIsNonlin_ || firstTime_)
	      // fetch solution at element nodes
	      //GetSolOfElement(elSol, connect_PDE);
	      sol_->GetElemSolutionAsMatrix(elSol, connect_PDE);
	      
	    
	    // ================================================================
	    //                             assemble matrices
	    // ================================================================

	    for(Integer actInteg=0; actInteg < integrators_[actDom]->size();
		actInteg++)
	      {
		IntegratorDescriptor * actDescriptor =
		  (*integrators_[actDom])[actInteg];


		// assemble only if nonlinear or first time
		if (reassembleMat_[actDescriptor->DestMat()] || firstTime_)
		  {
		    actDescriptor->GetIntegrator()->SetElemPtr(ptEl);
#ifdef USE_OLAS
		    FEMatrixType destMat = actDescriptor->DestMat();
#else
		    enum MatrixType destMat = actDescriptor->DestMat();
#endif
			    
		    // this matrix is nonlinear and, therefore, has to be reassembled next time
		    if (actDescriptor->IsNonLin())
		      {
			oneIntIsNonlin_ = TRUE;
			reassembleMat_[actDescriptor->DestMat()] = TRUE;
		
			actDescriptor->GetIntegrator()->SetActElemSol(elSol);
		      }

		    actDescriptor->GetIntegrator()->CalcElementMatrix(ptCoord, elemmat);
		    algsys_->SetElementMatrix(elemmat.GetDataPointer(), connect_PDE.GetPointer(), 
					      connect_PDE.GetSize(), destMat);
#ifdef DEBUG
		    (*debug) << "ElementMatrix of Element " << actEl << std::endl;
		    (*debug) << elemmat << std::endl;
#endif		    

		    if (actDescriptor->GetSecondaryMat() != NOTYPE)
		      {
			elemmat *= actDescriptor->GetSecMatFac();
			algsys_->SetElementMatrix(elemmat.GetDataPointer(), connect_PDE.GetPointer(), 
						  connect_PDE.GetSize(), actDescriptor->GetSecondaryMat()); 
		      }
		    
		  }		
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
	    if (oneIntIsNonlin_ || firstTime_)
	      // fetch solution at element nodes
	      //GetSolOfElement(elSol, connect_PDE);
	      sol_->GetElemSolutionAsMatrix(elSol, connect_PDE);
	      
	    
	    // ================================================================
	    //                             assemble matrices
	    // ================================================================

	    for(Integer actInteg=0; actInteg < surfintegrators_[actDom]->size(); actInteg++)
	      {
		IntegratorDescriptor * actDescriptor =
		  (*surfintegrators_[actDom])[actInteg];


		// assemble only if nonlinear or first time
		if (reassembleMat_[actDescriptor->DestMat()] || firstTime_)
		  {
		    actDescriptor->GetIntegrator()->SetElemPtr(ptEl);
#ifdef USE_OLAS
		    FEMatrixType destMat = actDescriptor->DestMat();
#else
		    enum MatrixType destMat = actDescriptor->DestMat();
#endif

		    // this matrix is nonlinear and, therefore, has to be reassembled next time
		    if (actDescriptor->IsNonLin())
		      {
			oneIntIsNonlin_ = TRUE;
			reassembleMat_[actDescriptor->DestMat()] = TRUE;
		
			// ?????????????????????????????????????????
			// was commented out before (18.3.04)
			actDescriptor->GetIntegrator()->SetActElemSol(elSol);
		      }

		    actDescriptor->GetIntegrator()->CalcElementMatrix(ptCoord, elemmat);
		    
		    algsys_->SetElementMatrix(elemmat.GetDataPointer(), connect_PDE.GetPointer(), 
					      connect_PDE.GetSize(), destMat); 

		    if (actDescriptor->GetSecondaryMat()  != NOTYPE )
		      {
			elemmat *= actDescriptor->GetSecMatFac();
			algsys_->SetElementMatrix(elemmat.GetDataPointer(), connect_PDE.GetPointer(), 
						  connect_PDE.GetSize(), actDescriptor->GetSecondaryMat()); 
		      }
		  }		
	      }
	  }
      }

     firstTime_ = FALSE;
  }




  // do the basic assembling stuff
  void Assemble::AssembleSrcRHS(const Integer level, const Double time)
  {
    ENTER_FCN( "Assemble:AssembleSrcRHS" );
    AssembleRHSNodalSources(level, time);
    AssembleRHSIntegralSources(level, time);
  }

  




  // do the basic assembling stuff
  void Assemble::AssembleRHSIntegralSources(const Integer level,
					    const Double time)
  {
    ENTER_FCN( "Assemble:AssembleRHSIntegralSources" );
 
     for (Integer actDom=0; actDom <  subdoms_.size(); actDom++)
      {	
	if (rhsSrcIntegrators_[actDom]->size())
	  {
	    std::vector<Elem*> elemssd;
	    ptgrid_->GetElemSD(elemssd, subdoms_[actDom], level);
	    
	    Double val_tfunc = 1.0;
	    if (ptTimeFunc_->GetmaxTimeFnc() > 0 )
		val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,fncname_rhs_[actDom]);

	    for (Integer actEl=0; actEl< elemssd.size(); actEl++)
	      {	       
		BaseFE * ptEl = elemssd[actEl]->ptElem;
		Vector<Integer> connecth = elemssd[actEl]->connect;
		
		Matrix<Double> ptCoord;
		GetElemCoords(connecth, ptCoord, level);
	    
		// map connect to PDE node numbers
		Vector<Integer> connect_PDE;
		Mesh2PDENode(connect_PDE, connecth, *mesh2PDENode_);
		
		for(Integer actRhsInt=0; actRhsInt < rhsSrcIntegrators_[actDom]->size(); actRhsInt++)
		  {
		    BaseIntDescriptor * actRhsID = (*rhsSrcIntegrators_[actDom])[actRhsInt];
		    
		    actRhsID->GetIntegrator()->SetElemPtr(ptEl);
		    
		    std::vector<Double> elemVec;
		    actRhsID->GetIntegrator()->CalcElemVector(ptCoord, elemVec);
		    
		    if (val_tfunc != 1.0)
		      elemVec *= val_tfunc;
		    
		    algsys_->SetElementRHS(&elemVec[0], connect_PDE.GetPointer(), connect_PDE.GetSize());
		  }
	      }
	  }
      }
  }



  // do the basic assembling stuff
  void Assemble::AssembleRHSNodalSources(const Integer level, const Double time)
  {
    ENTER_FCN( "Assemble:AssembleRHSNodalSources" );
    
    for (int actDom=0; actDom < loadDom_.size(); actDom++)
      {
	std::string doftype = loadDom_[actDom];

	Integer dof = 1;
	if (dofsPerNode_ != 1)
	  dof = GetBCDof( loadDof_[actDom] );	

	std::list<Integer> nodes;
	nodes = ptBCs_->GetNodesLevel(loadDom_[actDom], level);
	
	Double val = loadVals_[actDom];

	Double val_tfunc = 1.0;
	if (ptTimeFunc_->GetmaxTimeFnc() > 0 )
	    val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,fncname_loads_[actDom]);


	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
	  {
	    Integer node = *p;
	    
	    val = loadVals_[actDom] * val_tfunc;
	    algsys_->SetNodeRHS(val, (*mesh2PDENode_)[node-1], dof);	
	  }
      }
  }
  






  // do the basic assembling stuff
  void Assemble::AssembleNLRHS(const Integer level, const Double time)
  {
    ENTER_FCN( "Assemble:AssembleNLRHS" );

    Matrix<Double> elemmat;

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

	    sol_->GetElemSolutionAsMatrix(elSol, connect_PDE);
	      
	    
	    // ================================================================
	    //                             assemble RHS
	    // ================================================================

	    for(Integer actRhsInt=0; actRhsInt < rhsIntegrators_[actDom]->size(); actRhsInt++)
	      {
		BaseIntDescriptor * actRhsID = (*rhsIntegrators_[actDom])[actRhsInt];

		actRhsID->GetIntegrator()->SetElemPtr(ptEl);

		if (actRhsID->IsNonLin())
		  actRhsID->GetIntegrator()->SetActElemSol(elSol);

		
		std::vector<Double> elemVec;
		actRhsID->GetIntegrator()->CalcElemVector(ptCoord, elemVec);
		
		algsys_->SetElementRHS(&elemVec[0], connect_PDE.GetPointer(), connect_PDE.GetSize());
	      }
	  }
      }
  }


  Integer Assemble::
  GetBCDof(const std::string dofString)
  {
    ENTER_FCN( "MechPDE::GetBCDof" );

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
    ENTER_FCN( "Assemble::Mesh2PDENode" );
    PDENodes.Resize(MeshNodes.GetSize());
    
    for (Integer i=0; i<MeshNodes.GetSize(); i++) 
      PDENodes[i] = Mesh2PDENode[MeshNodes[i]-1];
  }





  void Assemble::InitMatrices()
  {
    ENTER_FCN( "Assemble::InitMatrices" );

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
 


  void Assemble::InitNonLinMatrices()
  {
    ENTER_FCN( "Assemble::InitNonLinMatrices" );

    // return, if matrices are not yet assembled
    if (!reassembleMat_.GetSize())
      {
	InitMatrices();
	return;
      }
    
    // Initialize matrices in order to get BCs correct
    algsys_->InitMatrix(SYSTEM);

    if (stiffnessMatrix_ && reassembleMat_[STIFFNESS])
      algsys_->InitMatrix(STIFFNESS);
    
    if (dampingMatrix_ && reassembleMat_[DAMPING])
      algsys_->InitMatrix(DAMPING);

    if (convectionMatrix_ && reassembleMat_[CONVECTION])
      algsys_->InitMatrix(CONVECTION);
    
    if (massMatrix_ && reassembleMat_[MASS])
      algsys_->InitMatrix(MASS); 
  }
 


 

  void Assemble::CreateMatrices()
  {
    ENTER_FCN( "Assemble::CreateMatrices" );
    const Integer numconstraints = 0;  // currently not handled
    
#ifdef USE_OLAS
    FEMatrixType matrixsystype[5];
    matrixsystype[0] = OLAS::NOTYPE;
    matrixsystype[1] = OLAS::NOTYPE;
    matrixsystype[2] = OLAS::NOTYPE;
    matrixsystype[3] = OLAS::NOTYPE;
    matrixsystype[4] = OLAS::NOTYPE;
#else
    Integer matrixsystype[5];
#endif
    
    matrixsystype[0] = SYSTEM;      // memory for the system matrix
    if (stiffnessMatrix_  == 1) matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
    if (dampingMatrix_    == 1) matrixsystype[2] = DAMPING;     // memory for the damping matrix
    if (convectionMatrix_ == 1) matrixsystype[3] = CONVECTION;  // memory for the convection matrix
    if (massMatrix_       == 1) matrixsystype[4] = MASS;        // memory for the mass matrix
    
    //put to algebraic system

#ifdef USE_OLAS
    olasParams_->SetValue( "FEMatrixType1", matrixsystype[0] );
    olasParams_->SetValue( "FEMatrixType2", matrixsystype[1] );
    olasParams_->SetValue( "FEMatrixType3", matrixsystype[2] );
    olasParams_->SetValue( "FEMatrixType4", matrixsystype[3] );
    olasParams_->SetValue( "FEMatrixType5", matrixsystype[4] );
    olasParams_->SetValue( "MatrixEntryType", entryType_ );
    olasParams_->SetValue( "MatrixStorageType", storageType_ );
    olasParams_->SetValue( "NumDof", dofsPerNode_ );
    olasParams_->SetValue( "NumDirichletBCs", numDirichletBCs_ );
    olasParams_->SetValue( "NumConstraints", numconstraints );
    olasParams_->SetValue( "AuxiliaryMatrix", FALSE);
    
    algsys_->CreateLinSys();
#else
    algsys_->CreateMatrix(matrixsystype, matrixType_, graphType_, 
			  dofsPerNode_,numDirichletBCs_, numconstraints);
#endif
    
  }
  

  void Assemble::SetGeneralParams(const std::string & pdename, 
				  const Integer dofsPerNode,
				  const Integer numPDENodes, 
				  const std::vector<std::string> subdoms,
				  const std::vector<std::string> surfdoms)
  {
    ENTER_FCN( "Assemble::SetGeneralParams" );

    pdename_     = pdename;
    dofsPerNode_ = dofsPerNode;
    numPDENodes_ = numPDENodes;
    subdoms_     = subdoms;
    surfdoms_    = surfdoms;



    // read load values =========================================

#ifndef XMLPARAMS
    conf->ifgetliststr("loads", loadDom_, pdename_);
    conf->ifgetliststr("loadDof", loadDof_, pdename_);


    if (dofsPerNode_ != 1)

      //check for load data
      if (loadDom_.size() != loadDof_.size())
	{
	  std::string errmsg = "Inconsistent definition of loads\n";
	  errmsg += "Dirichlet Boundary Conditions\n";
	  errmsg += " loadDom_.size() = " + Info->GenStr(loadDom_.size());
	  errmsg += "\n loadDof_.size() = " + Info->GenStr(loadDof_.size())
	    + '\n';
	  Info->Error( errmsg, __FILE__, __LINE__ );
	}
    
    loadVals_.resize(loadDom_.size());
    fncname_loads_.resize(loadDom_.size());

    for( int i = 0; i < loadDom_.size(); i++ )
      {
	conf->get2(loadDom_[i], loadVals_[i], fncname_loads_[i], pdename_,
		   "bc_conditions","loads");
      }

#else
    params->GetList( "name"    , loadDom_      , pdename_, "load" );
    params->GetList( "dof"     , loadDof_      , pdename_, "load" );
    params->GetList( "value"   , loadVals_     , pdename_, "load" );
    params->GetList( "dynamics", fncname_loads_, pdename_, "load" );

    // Check consistency
    if ( loadDom_.size() != loadDof_.size() ||
	 loadDom_.size() != loadVals_.size() )
      {
	std::string errmsg = "Loads: ";
	errmsg += "#name = " + Info->GenStr(loadDom_.size());
	errmsg += ", #dof = " + Info->GenStr(loadDof_.size());
	errmsg += ", #value = " + Info->GenStr(loadVals_.size());
	errmsg += ", #dynamics = " + fncname_loads_.size() + '\n';
	Info->Error( errmsg, __FILE__, __LINE__ );
      }

    // We need not have as many function/filenames as loads!
    for ( Integer k = fncname_loads_.size(); k < loadDom_.size(); k++ )
      {
	fncname_loads_.push_back( "none" );
      }
#endif

#ifdef DEBUG
    (*debug) << "Assemble::SetGeneralParams: We got " << loadDom_.size()
	     << " interfaces with loads" << std::endl;
    (*debug) << "Loads: #interfaces = " << loadDom_.size()
	     << ", #dof = " << loadDof_.size()
	     << ", #value = " << loadVals_.size()
	     << ", #dynamics = " << fncname_loads_.size() << std::endl;
    for ( unsigned int k = 0; k < loadDom_.size(); k++ )
      {
	(*debug) << "Loads: interface = " << loadDom_[k]
		 << ", dof = " << loadDof_[k]
		 << ", value = " << loadVals_[k];
	if ( k < fncname_loads_.size() )
	  {
	    (*debug) << ", dynamics = " << fncname_loads_[k] << std::endl;
	  }
	else
	  {
	    (*debug) << ", dynamics = " << std::endl;
	  }
      }
#endif

    // for every domain, we need an own integrator list ==========
    integrators_.resize(subdoms_.size());
    rhsSrcIntegrators_.resize(subdoms_.size());
    fncname_rhs_.resize(subdoms_.size());

    surfintegrators_.resize(surfdoms_.size());
    rhsIntegrators_.resize(subdoms_.size());
    for (int i=0; i<subdoms_.size();i++)
      {
	integrators_[i] = new std::vector<IntegratorDescriptor *>;
	rhsSrcIntegrators_[i] = new std::vector<BaseIntDescriptor *>;
	rhsIntegrators_[i] = new std::vector<BaseIntDescriptor *>;
      }

    for (int i=0; i<surfdoms_.size();i++)
      surfintegrators_[i] = new std::vector<IntegratorDescriptor *>;
  }
  

  //  void Assemble::SetupMatrixGraph(Integer numeq, Integer graphType)
  void Assemble::SetupMatrixGraph(Integer numeq)
  {
    ENTER_FCN( "Assemble::SetupMatrixGraph" );
    

  //initialize matrix graph
 
#ifdef USE_OLAS
  algsys_->CreateGraph(numeq, graphType_);
#else
  Integer estimated_maxneighbors=20;
  algsys_->CreateGraph(numeq, graphType_,estimated_maxneighbors);
#endif

  // set the graph - connectivity matrix
  BaseFE * ptElem; 
  Integer nsub, iel;

#ifdef USE_OLAS
  FEType fe_type;
#else
  Integer fe_type;
#endif
  
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
	  algsys_->SetElementPos(connecth.GetPointer(),connecth.GetSize(),fe_type);
	}
    }
  }



  Integer Assemble::GetNrBCDof(const std::string & dofStartString)
  {
    ENTER_FCN( "Analysis::GetNrBCDof" );
    
    Integer nrActDof = 0;
    
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
  


  // define integrators
  void Assemble::AddRhsIntegrator(BaseForm * integrator,
				  const std::string & subDomName, 
				  const Integer nonLin)
  {
    ENTER_FCN( "Assemble::AddRhsIntegrator" );
    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsIntegrators_[SubDomIndex(subDomName)]->push_back(actRhsID);
  }


  void Assemble::AddRhsSrcIntegrator(BaseForm * integrator,
				     const std::string & subDomName, 
				     const std::string fncname,
				     const Integer nonLin)
  {
    ENTER_FCN( "Assemble::AddRhsSrcIntegrator" );
    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsSrcIntegrators_[SubDomIndex(subDomName)]->push_back(actRhsID);
    fncname_rhs_[SubDomIndex(subDomName)] = fncname;
  }


  // ==========================================================
  // STATIC ANALYSIS
  // ==========================================================


  StaticAssemble::StaticAssemble(BaseSystem * algsys, Grid * agrid)
    :Assemble(algsys, agrid)
  {
    ENTER_FCN( "StaticAssemble::StaticAssemble" );
    graphType_    = NODEGRAPH; 
  }
  

#ifdef USE_OLAS 
  /// define integrators
  void StaticAssemble::AddIntegrator(BaseForm * integrator, const std::string & subDomName,
					const FEMatrixType destinationMatrix, const Integer nonLin)
  {
#ifdef TRACE
    (*trace) << "entering StaticAssemble::AddIntegrator " << std::endl;
#endif

#ifdef USE_OLAS
    FEMatrixType actMatType = destinationMatrix;
#else
    MatrixType actMatType = destinationMatrix;
#endif

    
    if (actMatType == STIFFNESS)
      actMatType = SYSTEM;

    if (actMatType !=  SYSTEM)
      return;

    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator, actMatType, nonLin);
    integrators_[SubDomIndex(subDomName)]->push_back(actID);

  }

#else
  // define integrators
  void StaticAssemble::AddIntegrator(BaseForm * integrator,
				     const std::string & subDomName,
				     const enum MatrixType destinationMatrix,
				     const Integer nonLin)
  {
    ENTER_FCN( "StaticAssemble::AddIntegrator" );

    MatrixType actMatType = destinationMatrix;
    
    if (actMatType == STIFFNESS)
      actMatType = SYSTEM;

    if (actMatType !=  SYSTEM)
      return;

    IntegratorDescriptor * actID =
      new IntegratorDescriptor(integrator, actMatType, nonLin);
    integrators_[SubDomIndex(subDomName)]->push_back(actID);
  }


#endif

#ifdef USE_OLAS 
  /// define integrators
  void StaticAssemble::AddSurfIntegrator(BaseForm * integrator, const std::string & subDomName,
					const FEMatrixType destinationMatrix, const Integer nonLin)
  {
#ifdef TRACE
    (*trace) << "entering StaticAssemble::AddSurfIntegrator " << std::endl;
#endif
   
#ifdef USE_OLAS
    FEMatrixType actMatType = destinationMatrix;
#else
    MatrixType actMatType = destinationMatrix;
#endif
    
    if (actMatType == STIFFNESS)
      actMatType = SYSTEM;

    if (actMatType !=  SYSTEM)
      return;

    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator, actMatType, nonLin);
    surfintegrators_[SurfDomIndex(subDomName)]->push_back(actID);
  }
#else
  /// define integrators
  void
  StaticAssemble::AddSurfIntegrator(BaseForm * integrator,
				    const std::string & subDomName,
				    const enum MatrixType destinationMatrix,
				    const Integer nonLin)
  {
    ENTER_FCN( "StaticAssemble::AddSurfIntegrator" );

    MatrixType actMatType = destinationMatrix;
    
    if (actMatType == STIFFNESS)
      actMatType = SYSTEM;

    if (actMatType !=  SYSTEM)
      return;

    IntegratorDescriptor * actID =
      new IntegratorDescriptor(integrator, actMatType, nonLin);

    surfintegrators_[SurfDomIndex(subDomName)]->push_back(actID);
  }
#endif

  // define integrators
  void StaticAssemble::AddIntegrator(IntegratorDescriptor * actID,
				     const std::string & subDomName)
  {
    ENTER_FCN( "StaticAssemble::AddIntegrator" );

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
    ENTER_FCN( "TransientAssemble::TransientAssemble" );
    graphType_    = NODEGRAPH; 
    stiffnessMatrix_  = TRUE;
    massMatrix_       = TRUE;
  }


#ifdef USE_OLAS 
  /// define integrators
  void TransientAssemble::AddIntegrator(BaseForm * integrator, const std::string & subDomName,
					const FEMatrixType destinationMatrix, const Integer nonLin)
  {
#ifdef TRACE
    (*trace) << "entering TransientAssemble::AddIntegrator " << std::endl;
#endif
    if (destinationMatrix == SYSTEM)
      Info->Error("In transient assembling, no SYSTEM matrix may be defined directly", __FILE__, __LINE__);

    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator, destinationMatrix, nonLin);
    integrators_[SubDomIndex(subDomName)]->push_back(actID);
  }
#else
  /// define integrators
  void TransientAssemble::AddIntegrator(BaseForm * integrator, const std::string & subDomName,
					const enum MatrixType destinationMatrix, const Integer nonLin)
  {
#ifdef TRACE
    (*trace) << "entering TransientAssemble::AddIntegrator " << std::endl;
#endif
    if (destinationMatrix == SYSTEM)
      Info->Error("In transient assembling, no SYSTEM matrix may be defined directly", __FILE__, __LINE__);

    IntegratorDescriptor * actID =
      new IntegratorDescriptor(integrator, destinationMatrix, nonLin);
    integrators_[SubDomIndex(subDomName)]->push_back(actID);
  }
#endif


#ifdef USE_OLAS
  /// define integrators
  void TransientAssemble::AddSurfIntegrator(BaseForm * integrator, const std::string & subDomName,
					const FEMatrixType destinationMatrix, const Integer nonLin)
  {
#ifdef TRACE
    (*trace) << "entering TransientAssemble::AddSurfIntegrator " << std::endl;
#endif
    if (destinationMatrix == SYSTEM)
      Info->Error("In transient assembling, no SYSTEM matrix may be defined directly", __FILE__, __LINE__);

    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator, destinationMatrix, nonLin);
    surfintegrators_[SurfDomIndex(subDomName)]->push_back(actID);
  }
#else


  // define integrators
  void
  TransientAssemble::AddSurfIntegrator(BaseForm * integrator,
				       const std::string & subDomName,
				       const enum MatrixType destinationMatrix,
				       const Integer nonLin)
  {
    ENTER_FCN( "TransientAssemble::AddSurfIntegrator" );

    if (destinationMatrix == SYSTEM)
      Info->Error("In transient assembling, no SYSTEM matrix may be defined directly", __FILE__, __LINE__);
    
    IntegratorDescriptor * actID =
      new IntegratorDescriptor(integrator, destinationMatrix, nonLin);
    surfintegrators_[SurfDomIndex(subDomName)]->push_back(actID);
  }
#endif




  // define integrators
  void TransientAssemble::AddIntegrator(IntegratorDescriptor * actID,
					const std::string & subDomName)
  {
    ENTER_FCN( "TransientAssemble::AddIntegrator" );

    if (actID->DestMat() == SYSTEM)
      Info->Error("In transient assembling, no SYSTEM matrix may be defined directly", __FILE__, __LINE__);

    integrators_[SubDomIndex(subDomName)]->push_back(actID);
  }



  // ==========================================================
  // RHS INTEGRATOR DESCRIPTOR
  // ==========================================================

  BaseIntDescriptor::BaseIntDescriptor() : integrator(NULL), nonLin(FALSE)
  {
    ENTER_FCN( "BaseIntDescriptor::BaseIntDescriptor" );
  }

  BaseIntDescriptor::~BaseIntDescriptor()
  {
    ENTER_FCN( "BaseIntDescriptor::~BaseIntDescriptor" );
    if (integrator)
      delete integrator;
  }
  

  // define integrators
  BaseIntDescriptor::BaseIntDescriptor(BaseForm * aIntegrator,
				       const Boolean aNonLin)
    : integrator(aIntegrator), nonLin(aNonLin)
  {
    ENTER_FCN( "BaseIntDescriptor::BaseIntDescriptor " );
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
    ENTER_FCN( "IntegratorDescriptor::IntegratorDescriptor" );
  }


#ifdef USE_OLAS
  /// define integrators
  IntegratorDescriptor::IntegratorDescriptor(BaseForm * aIntegrator, 
					     const enum FEMatrixType aDestMat, const Boolean aNonLin)
    :BaseIntDescriptor(aIntegrator, aNonLin),
     destinationMatrix(aDestMat),
     secondaryMatrix(NOTYPE),
     secMatFac(0.0)
{
#ifdef TRACE
  (*trace) << "entering IntegratorDescriptor::IntegratorDescriptor" << std::endl;
#endif
  }
  
#else
  /// define integrators
  IntegratorDescriptor::IntegratorDescriptor(BaseForm * aIntegrator, 
					     const enum MatrixType aDestMat,
					     const Boolean aNonLin)
    :BaseIntDescriptor(aIntegrator, aNonLin),
     destinationMatrix(aDestMat),
     secondaryMatrix(NOTYPE),
     secMatFac(0.0)
  {
    ENTER_FCN( "IntegratorDescriptor::IntegratorDescriptor" );
  }
  
#endif


  IntegratorDescriptor::~IntegratorDescriptor()
  {
    ENTER_FCN( "IntegratorDescriptor::~IntegratorDescriptor" );
    if (integrator)
      delete integrator;
  }


} // end namespace CoupledField
