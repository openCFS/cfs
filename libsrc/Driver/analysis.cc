#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include <DataInOut/conffile.hh>
#include <DataInOut/WriteInfo.hh>
#include <multigrid.hh>

#include "analysis.hh" 

namespace CoupledField
{



  Analysis::Analysis(Grid * aptgrid)
    :ptgrid_(aptgrid),
     stiffnessMatrix_(FALSE),
     dampingMatrix_(FALSE),
     massMatrix_(FALSE),
     convectionMatrix_(FALSE),
     actlevel_(0),
     integrators_(0)

  {
#ifdef TRACE
    (*trace) << "entering Analysis::Analysis " << std::endl;
#endif
  }





  Analysis::~Analysis()
  {
#ifdef TRACE
    (*trace) << "entering Analysis::~Analysis " << std::endl;
#endif
    
    for (int i=0; i<integrators_.size();i++)
      for (int j=0; j<integrators_[i]->size(); j++)
	delete (*integrators_[i])[j];
  }




   

 
  /// define integrators
  void Analysis::AddIntegrator(BaseForm * integrator, const std::string & subDomName,
			       const enum MatrixType destinationMatrix, const Integer nonLin)
  {
#ifdef TRACE
    (*trace) << "entering Analysis::AddIntegrator " << std::endl;
#endif

    struct integratorDescriptor * actID = new struct integratorDescriptor;
    
    actID->integrator = integrator;
    actID->destinationMatrix = destinationMatrix;
    actID->nonLin = nonLin;

    integrators_[SubDomIndex(subDomName)]->push_back(actID);
  }
  







  Integer Analysis::SubDomIndex(const std::string & subDomName)
  {
#ifdef TRACE
    (*trace) << "entering Analysis::SubDomIndex " << std::endl;
#endif



    for (int i=0; i < subdoms_.size();i++)
      if (subDomName == subdoms_[i])
	return i;
  
  
    std::string errOut;
    errOut = "SubDomain " + subDomName + " not defined!";
    Info->Error(errOut, __FILE__, __LINE__);
  }
  






  void Analysis::GetElemCoords(const Vector<Integer> connect, 
			       Matrix<Double> &coordMat, const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering Analysis:GetElemCoords:" << std::endl;
#endif
    
    ptgrid_->GetCoordNodesElemMat(connect, coordMat, level);
  
    //    if (deltCoords_.size() != 0)
    //      {
    //     for (Integer i=0; i<coordMat.size_row(); i++)
    //       for (Integer j=0; j<coordMat.size_col(); j++) 
    // 	coordMat(i,j) += deltCoords_[i][Mesh2PDENode_ [connect[j] - 1] - 1];
    //      }
  }







  
  
  /// define integrators
  void Analysis::SetupMatrices(const Integer level)
  {
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


	    for(int actInteg=0; actInteg < integrators_[actDom]->size(); actInteg++)
	      {
		struct integratorDescriptor * actDescriptor = (*integrators_[actDom])[actInteg];
		enum MatrixType destMat = actDescriptor->destinationMatrix;
		
		actDescriptor->integrator->CalcElementMatrix(ptCoord, elemmat);
		algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), destMat);
	      }
	  }
      }
  }





  void Analysis::Mesh2PDENode(Vector<Integer> & PDENodes, 
			      const Vector<Integer> & MeshNodes,
			      const std::vector<Integer> & Mesh2PDENode)
  {
#ifdef TRACE
    (*trace) << "entering Analysis::Mesh2PDENode " << std::endl;
#endif
    
    PDENodes.Resize(MeshNodes.size());
    
    for (Integer i=0; i<MeshNodes.size(); i++) 
      PDENodes[i] = Mesh2PDENode[MeshNodes[i]-1];
  }  



  void Analysis::CreateMatrices_Solver()
  {
#ifdef TRACE
    (*trace) << "entering  Analysis::CreateMatrices_Solver" << std::endl;
#endif
    
    Integer matrixsystype[5];
  
    // system matrix is always needed
    matrixsystype[0] = SYSTEM;                                  // memory for the system matrix
    if (stiffnessMatrix_  == 1) matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
    if (dampingMatrix_    == 1) matrixsystype[2] = DAMPING;     // memory for the damping matrix
    if (convectionMatrix_ == 1) matrixsystype[3] = CONVECTION;  // memory for the convection matrix
    if (massMatrix_       == 1) matrixsystype[4] = MASS;        // memory for the mass matrix

    Integer numConstraints = 0;  // currently not handled

    //put to algebraic system
    algsys_->CreateMatrix(matrixsystype, matrixType_, graphType_, dofsPerNode_, numDirichletBCs_,
			  numConstraints);

    //create solver and preconditioner
    algsys_->CreateSolver();
    algsys_->CreatePrecond(matrixType_);

    //now reset AlgebraicSystem 
    algsys_->InitRHS();
    algsys_->InitSol();

    InitMatrices();  
  }




  void Analysis::InitMatrices()
  {
    // Initialize matrices in order to get BCs correct
    if (stiffnessMatrix_)
      algsys_->InitMatrix(STIFFNESS);
    
    if (dampingMatrix_)
      algsys_->InitMatrix(DAMPING);

    if (convectionMatrix_)
      algsys_->InitMatrix(CONVECTION);
    
    if (massMatrix_)
      algsys_->InitMatrix(MASS); 
  }
  


  void Analysis::SetGeneralParams(const std::string & pdename, 
				  const Integer dofsPerNode,
				  const Integer numPDENodes, 
				  const std::vector<std::string> subdoms)
  {
    pdename_     = pdename;
    dofsPerNode_ = dofsPerNode;
    numPDENodes_ = numPDENodes;
    subdoms_     = subdoms;

    // for every domain, we need an own integrator list
    integrators_.resize(subdoms_.size());
    for (int i=0; i<subdoms_.size();i++)
      //      integrators_[i] = new Vector<integratorDescriptor *>;
      integrators_[i] = new std::vector<integratorDescriptor *>;
    
  }
  


  void Analysis::SetSolverParameters()
  {
#ifdef TRACE
    (*trace) << "entering  BasePDE::SetSolverParameters" << std::endl;
#endif


    //solver parameters
    Integer maxnumiter  = 100;
    Integer solvertype  = RealDirect;
    Integer precondtype = ID;
    Integer numeqcoarse = 200;
    Double  eps         = 1e-8;
    Double  dampiter    = 1.0;
    Double  coarsealpha = 0.1;

    //if values are defined in conf-file, take these
    conf->ifget("eps", eps, pdename_); // relative accuracy in the precond. energy
    conf->ifget("dampiter", dampiter, pdename_); // damping parameter for Jacobi, SSOR
    conf->ifget("maxnumit", maxnumiter, pdename_); // maximal number of iterations
    conf->ifget("solvertype", solvertype, pdename_); // solver
    conf->ifget("precondtype", precondtype, pdename_); //preconditioner
    conf->ifget("numeqcoarse", numeqcoarse, pdename_); // number of equation for coarsing
    conf->ifget("coarsealpha", coarsealpha, pdename_); // coarsing parameter for AMG
  
    if (solvertype==RealDirect && precondtype!=ID) 
      precondtype=ID;
  
    //communicate with algebraic system
    algsys_->CreateParameter();
    algsys_->SetAccuracy(eps);
    algsys_->SetMaxNumIter(maxnumiter);
    algsys_->SetPrecond(precondtype);
    algsys_->SetSolver(solvertype);
    algsys_->SetDampIter(dampiter);
    algsys_->SetCoarseSystem(numeqcoarse);
    algsys_->SetAlpha(coarsealpha);

  } 



  void Analysis::SetAlgSys()
  {
#ifdef TRACE
    (*trace) << "entering  Analysis::SetAlgSys" << std::endl;
#endif

    //allocate according algebraic system
    algsys_ = new StandardSystem();

    //set solver parameters  
    SetSolverParameters();

    //ask the PDE matrix settings
    MatrixSettings();

    //set the graph type used for the system matrices
    SetupMatrixGraph(numPDENodes_, graphType_);

    //allocate the necessary matrices as well as solver and preconditioner
    CreateMatrices_Solver();
  }






void Analysis::SetupMatrixGraph(Integer numeq, Integer graphtype)
{
#ifdef TRACE
  (*trace) << "entering Analysis::SetupMatrixGraph" << std::endl;
#endif

  //initialize matrix graph
  algsys_->CreateGraph(numeq,graphtype);

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


    void Analysis::DeleteAlgSys(){if (algsys_) delete algsys_;};

  // ==========================================================
  // STATIC ANALYSIS
  // ==========================================================

  void StaticAnalysis::MatrixSettings()
  {
    matrixType_   = RBLOCK;
    graphType_    = NODEGRAPH; 
  }


  // "time stepping" for solver
  void StaticAnalysis::SolveStep()
  {
  };
  

} // end namespace CoupledField
