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
  
  Assemble::Assemble(BaseSystem * algsys, Grid * aptgrid)
    :algsys_(algsys),
     ptgrid_(aptgrid),
     stiffnessMatrix_(FALSE),
     dampingMatrix_(FALSE),
     massMatrix_(FALSE),
     convectionMatrix_(FALSE),
     actlevel_(0),
     integrators_(0)

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
  }




   

 
  /// define integrators
  void Assemble::AddIntegrator(BaseForm * integrator, const std::string & subDomName,
			       const enum MatrixType destinationMatrix, const Integer nonLin)
  {
#ifdef TRACE
    (*trace) << "entering Assemble::AddIntegrator " << std::endl;
#endif

    struct integratorDescriptor * actID = new struct integratorDescriptor;
    
    actID->integrator = integrator;
    actID->destinationMatrix = destinationMatrix;
    actID->nonLin = nonLin;

    integrators_[SubDomIndex(subDomName)]->push_back(actID);
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
    static int setupMatricesCount = 1;
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
		struct integratorDescriptor * actDescriptor =
		  (*integrators_[actDom])[actInteg];


		// assemble only if nonlinear needs or first time
		if (actDescriptor->nonLin || setupMatricesCount == 1)
		  {
		    actDescriptor->integrator->SetElemPtr(ptEl);

		    enum MatrixType destMat = actDescriptor->destinationMatrix;
		
		    actDescriptor->integrator->CalcElementMatrix(ptCoord, elemmat);
		    algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), 
					      connect_PDE.size(), destMat); 
		  }		
	      }
	  }
      }
  }




  
  // do the basic assembling stuff
  void Assemble::AssembleRHS(const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering Assemble:AssembleRHS" << std::endl;
#endif
    
    for (int actDom=0; actDom < loadDoms_.size(); actDom++)
    {
	std::string doftype = loadDoms_[actDom];
	Integer dof = GetNrBCDof (doftype.substr(0,2));

	std::list<Integer> nodes;
	nodes = ptBCs_->GetNodesLevel(loadDoms_[actDom], level);

	Double val = val_loads_[actDom];

	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
	  {
	    Integer node = *p;

	    val = val_loads_[actDom];
	    algsys_->SetNodeRHS(val, (*mesh2PDENode_)[node-1], dof);	
	  }


    }
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
// //  //Initialize matrices in order to get BCs correct
//   std::vector<Integer> matrixsystype(5,0);    
//   if (systemMatrix_     == 1) matrixsystype[0] = SYSTEM;      // memory for the system matrix
//   if (stiffnessMatrix_  == 1) matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
//   if (dampingMatrix_    == 1) matrixsystype[2] = DAMPING;     // memory for the damping matrix
//   if (convectionMatrix_ == 1) matrixsystype[3] = CONVECTION;  // memory for the convection matrix
//   if (massMatrix_       == 1) matrixsystype[4] = MASS;        // memory for the mass matrix
  

//   for (Integer i=0;i<5;i++)
//     if (matrixsystype[i] !=0)
//       algsys_->InitMatrix(i+1);

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
				  const std::vector<std::string> subdoms)
  {
    pdename_     = pdename;
    dofsPerNode_ = dofsPerNode;
    numPDENodes_ = numPDENodes;
    subdoms_     = subdoms;

    // for every domain, we need an own integrator list
    integrators_.resize(subdoms_.size());
    for (int i=0; i<subdoms_.size();i++)
      integrators_[i] = new std::vector<integratorDescriptor *>;
      //      integrators_[i] = new Vector<integratorDescriptor *>;
    
  }
  












  //  void Assemble::SetupMatrixGraph(Integer numeq, Integer graphType)
  void Assemble::SetupMatrixGraph(Integer numeq)
  {
#ifdef TRACE
    (*trace) << "entering Assemble::SetupMatrixGraph" << std::endl;
#endif
    
  //initialize matrix graph
  algsys_->CreateGraph(numeq, graphType_);

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



  // ==========================================================
  // STATIC ANALYSIS
  // ==========================================================


  StaticAssemble::StaticAssemble(BaseSystem * algsys, Grid * agrid)
    :Assemble(algsys, agrid)
  {
    //    matrixType_   = RBLOCK;
    graphType_    = NODEGRAPH; 
  }
  



  // ==========================================================
  // TRANSIENT ANALYSIS
  // ==========================================================


  TransientAssemble::TransientAssemble(BaseSystem * algsys, Grid * agrid)
    :Assemble(algsys, agrid)
  {
    //    matrixType_   = RBLOCK;
    graphType_    = NODEGRAPH; 

    stiffnessMatrix_  = TRUE;
    massMatrix_       = TRUE;
  }
  


} // end namespace CoupledField
