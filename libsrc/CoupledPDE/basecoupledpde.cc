#include <fstream>
#include <iostream>
#include <string>

#include "basecoupledpde.hh"
#include <CoupledPDE/pdecoupling.hh>
#include <DataInOut/conffile.hh>
 
namespace CoupledField
{

BaseCoupledPDE::BaseCoupledPDE(std::vector<BasePDE*> & PDEs,
			       std::vector<PDECoupling*> & Couplings,
			       Grid *aptgrid, 
			       BCs *aptBCs, 
			       FileType *aInFile, 
			       WriteResults * aOutFile)
{
#ifdef TRACE
  (*trace) << "entering BaseCoupledPDE::BaseCoupledPDE" << std::endl;
#endif

  InFile_     = aInFile;
  OutFile_    = aOutFile;
  ptgrid_     = aptgrid;
  ptBCs_      = aptBCs;
  PDEs_       = PDEs;
  Couplings_  = Couplings;

  actlevel_ = 0;
  
  // get analysis type
  std::string analysis;
  conf->get("analysis", analysis);


  if (analysis=="static") 
    analysistype_ = STATIC;
  else if (analysis=="transient")
    analysistype_ = TRANSIENT;
  else
    Error("Analysis Type not supported",__FILE__,__LINE__);

  // set standard solver parameters
  eps_ = 1e-6;
  maxnumiter_ = 10;


  // FOR DEBUGGING ONLY
  
}



void BaseCoupledPDE::SetSolverParameters()
{
#ifdef TRACE
  (*trace) << "entering  BaseCoupledPDE::SetSolverParameters" << std::endl;
#endif

  //if values are defined in conf-file, take these
  conf->ifget("eps",eps_,pdename_); // relative accuracy in the precond. energy
  conf->ifget("maxnumit",maxnumiter_,pdename_); // maximal number of iterations
  
} 



void BaseCoupledPDE::Mesh2PDENode(Vector<Integer> & PDENodes, 
			   const Vector<Integer> & MeshNodes,
			   const std::vector<Integer> & Mesh2PDENode)
{
#ifdef TRACE
  (*trace) << "entering BaseCoupledPDE::Mesh2PDENode " << std::endl;
#endif

  PDENodes.Resize(MeshNodes.size());
  
  for (Integer i=0; i<MeshNodes.size(); i++)
    PDENodes[i] = Mesh2PDENode[MeshNodes[i]-1];

}


void BaseCoupledPDE::PDE2MeshNode(Vector<Integer> & MeshNodes, 
			   const Vector<Integer> & PDENodes,
			   const std::vector<Integer> & PDE2MeshNode)
{
#ifdef TRACE
  (*trace) << "entering BaseCoupledPDE::PDE2MeshNode " << std::endl;
#endif

  MeshNodes.Resize(PDENodes.size());

  for (Integer i=0; i<PDENodes.size(); i++)
    MeshNodes[i] = PDE2MeshNode[PDENodes[i]-1];

}

void BaseCoupledPDE::AssignPDENodeNumbers(std::vector<Integer> & Mesh2PDENode,
			  std::vector<Integer> & PDE2MeshNode,
			  const std::vector<std::string> & subdoms)
{
#ifdef TRACE
  (*trace) << "entering BaseCoupledPDE:AssignPDENodeNumbers:" << std::endl;
#endif

  // Initialize Mesh2PDENode and PDE2MeshNode
  Mesh2PDENode.resize(ptgrid_->GetMaxnumnodes(actlevel_),-1);
  std::vector<Elem*> SD;
  Integer NodeCounter = 1;
  
  // Iterate over Subdomains
  for (Integer i=0; i<subdoms.size(); i++)
    {
      ptgrid_->GetElemSD(SD,subdoms[i],actlevel_);
      // Iterate over all elements in subdomain
      for (Integer j=0; j<SD.size(); j++)
	{
	  // Iterate over all element nodes
	  for (Integer NumNodes=0; NumNodes<SD[j]->connect.size(); NumNodes++)
	    {
	      // Check if node was already assigned
	      if (Mesh2PDENode[SD[j]->connect[NumNodes] - 1] == -1)
		{
		  Mesh2PDENode[SD[j]->connect[NumNodes] - 1] = NodeCounter;
		  PDE2MeshNode.push_back(SD[j]->connect[NumNodes]);
		  NodeCounter++;
		}
	    }
	}
    }
}

  
void BaseCoupledPDE::AssignPDENodeNumbers(std::vector<Integer> & Mesh2PDENode,
			  std::vector<Integer> & PDE2MeshNode,
			  const std::vector<Elem*> & Elements)
{
#ifdef TRACE
  (*trace) << "entering BaseCoupledPDE:AssignPDENodeNumbers:" << std::endl;
#endif

  // Initialize Mesh2PDENode_ and PDE2MeshNode
  Mesh2PDENode.resize(ptgrid_->GetMaxnumnodes(actlevel_),-1);
  Integer NodeCounter = 1;
  
  // Iterate over all elements 
  for (Integer j=0; j<Elements.size(); j++)
    {
      // Iterate over all element nodes
      for (Integer NumNodes=0; NumNodes<Elements[j]->connect.size(); NumNodes++)
	{
	  // Check if node was already assigned
	  if (Mesh2PDENode[Elements[j]->connect[NumNodes] - 1] == -1)
	    {
	      Mesh2PDENode[Elements[j]->connect[NumNodes] - 1] = NodeCounter;
	      PDE2MeshNode.push_back(Elements[j]->connect[NumNodes]);
	      NodeCounter++;
	    }
	}
    }
}




void BaseCoupledPDE::TransformNodeSolution(Vector<Double> & MeshSol, 
				    const Vector<Double> & PDESol, 
				    const std::vector<Integer> & PDE2MeshNode,
				    const Integer nrDofs)
{
#ifdef TRACE
  (*trace) << "entering BaseCoupledPDE::TransformNodeSolution" << std::endl;
#endif

  Integer k, node, idx;
  MeshSol.Resize(ptgrid_->GetMaxnumnodes(actlevel_) * nrDofs);

  k=0;
  
  // Loop over all PDE nodes
  for (Integer i=0; i<PDE2MeshNode.size(); i++)
    {
      node = PDE2MeshNode[i];
      idx  = nrDofs *(node-1);
      for (Integer j=0; j<nrDofs; j++)
	{
	  MeshSol[idx+j] = PDESol[k];
	  k++;
	}
    }
}


void BaseCoupledPDE::TransformElemSolution(Vector<Double> & MeshSol, 
				    const Vector<Double> & PDESol, 
				    const std::vector<Elem*> & Elems,
				    const Integer nrDofs)
{
#ifdef TRACE
  (*trace) << "entering BaseCoupledPDE::TransformElemSolution" << std::endl;
#endif

  Integer ElemNum, k;
  
  MeshSol.Resize(ptgrid_->GetMaxnumElem(actlevel_) * nrDofs);
  
  k=0;
  
  // loop over all elements
  for (Integer i=0; i<Elems.size(); i++)
    {
      ElemNum = Elems[i]->ElemNum;
      
      for( Integer j=0; j<nrDofs; j++)
	{
	  MeshSol[ (ElemNum-1) * nrDofs+j] = PDESol[k++];
	}
    }
	      
}
void BaseCoupledPDE::TransformElemSolution(Vector<Double> & MeshSol, 
				    const Vector<Double> & PDESol, 
				    const std::vector<std::string> & SD,
				    const Integer nrDofs)
{
#ifdef TRACE
  (*trace) << "entering BaseCoupledPDE::TransformElemSolution" << std::endl;
#endif

  Integer ElemNum, k;
  std::vector<Elem*> Elems;

  MeshSol.Resize(ptgrid_->GetMaxnumElem(actlevel_) * nrDofs);

  k=0;

  // loop over all SubDomains
  for (Integer isd=0; isd<SD.size(); isd++)
  {
    ptgrid_->GetElemSD(Elems, SD[isd], actlevel_);
    //std::cerr << "size of SD " << SD[isd] << ":" << Elems.size() << std::endl;
    
    // loop over all elements
    for (Integer i=0; i<Elems.size(); i++)
      {
	ElemNum = Elems[i]->ElemNum - 1;
	
	for( Integer j=0; j<nrDofs; j++)
	  MeshSol[ElemNum*nrDofs+j] = PDESol[k++];
	
      }
  }
}

BaseCoupledPDE::~BaseCoupledPDE()
{
#ifdef TRACE
  (*trace) << " entering BaseCoupledPDE::~BaseCoupledPDE() " << std::endl;
#endif

}


} // end of namespace
