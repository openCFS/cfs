#include "pdecoupling.hh"

#include <Domain/grid.hh>
#include <Domain/bcs.hh>
#include <Domain/elem.hh>
#include <list>

#ifndef NEWBASEPDE
#include <PDE/basepde.hh>
#else
#include <PDE/newBasePDE.hh>
#endif //#ifndef NEWBASEPDE

namespace CoupledField
{

PDECoupling::PDECoupling(Grid * aptgrid, BCs * aptBCs)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::PDECoupling" << std::endl;
#endif 

  ptGrid_ = aptgrid;
  ptBCs_ = aptBCs;
  
  defaultEpsilon = 1e-5;
  defaultNormType = L2REL;

}

PDECoupling::~PDECoupling()
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::~PDECoupling" << std::endl;
#endif 
  
  for (Integer i=0; i<outputInterfaces_.size(); i++)
    if (outputInterfaces_[i]) delete outputInterfaces_[i];

   for (Integer i=0; i<inputInterfaces_.size(); i++)
    if (inputInterfaces_[i]) delete inputInterfaces_[i];

}

void PDECoupling::RegisterInput(CouplingInputType InType, std::string Quantity)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::RegisterInput" << std::endl;
#endif 
  
  inputTypes_.push_back(InType);
  inputQuantities_.push_back(Quantity);  
  inputInterfaces_.push_back(0);
}

void PDECoupling::AddInput(std::string quantity, 
			   std::string region, 
			   CouplingRegionType regionType,
			   Integer level,
			   std::vector<PDECoupling*> & couplings)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::AddInput" << std::endl;
#endif 

  //std::cerr << "AddInput called with Quantity \"" << Quantity << "\", region \"" << region << "\" regiontype \"" << RegionType << "\" and a vector of Couplings..." << std::endl;

  CouplingInterface *myInterface = 0; 
  Integer myNum = -1;

  // search matching  Quantity
  for (Integer i=0; i<inputQuantities_.size(); i++)
    {
      if (inputQuantities_[i] == quantity)
	myNum = i; 
    }

   if (myNum == -1)
    {
      std::string ErrMsg = "Qantity \'" + quantity + "\' not registered for PDE \'" + myPDE_->GetName()+ "\'";
      Error(ErrMsg.c_str(),__FILE__,__LINE__);
    }

 
  // Add coupling input as output to pde, which can calculate specified quantity
  CouplingOutputType myOutputType = Input2OutputType(inputTypes_[myNum]);

  Integer i=0;
  while (myInterface == 0 && i<couplings.size())
    {
      myInterface = couplings[i++]->AddOutput(myOutputType, quantity, region, regionType, level);
    }
  
  // If no pde has the specified quantity as output
  if (myInterface == NULL)
    {
      std::string ErrMsg = "Qantity \'" + quantity + "\' can not be calculated with current set of PDEs";
      Error(ErrMsg.c_str(),__FILE__,__LINE__);
    }
  
  
  // set normtype and epsilon
  myInterface->epsilon = defaultEpsilon;
  conf->ifget(inputQuantities_[myNum],myInterface->epsilon,"coupling", "tolerance"); 
  
  std::string normtype;
  myInterface->normtype = defaultNormType;
  
  if (conf->ifget(inputQuantities_[myNum],normtype,"coupling","normtype") == TRUE)
    {

      if (normtype == "L2rel")
	myInterface->normtype = L2REL;
      else if (normtype == "L2abs")
	myInterface->normtype = L2ABS;
      else 
	{
	  std::string errMsg = "Normtype \'" + normtype + "\' not known!";
	  Error(errMsg.c_str(),__FILE__,__LINE__);
	} 
    }
  
  inputInterfaces_[myNum] = myInterface;
          
  
}

PDECoupling::CouplingInterface* PDECoupling::AddOutput(CouplingOutputType outputType, 
						       std::string quantity, 
						       std::string region,
						       CouplingRegionType regionType,
						       Integer level)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::AddOutput" << std::endl;
#endif 

  // search if output exists already
  for (Integer i=0; i<outputTypes_.size(); i++)
    if (outputTypes_[i] == outputType 
	&& outputQuantities_[i] == quantity
	&& outputInterfaces_[i]->region == region)
      return outputInterfaces_[i];
  
  
  // if not find out if PDE can calculate desired output
  if (myPDE_->HasOutput(quantity) == FALSE)
    return 0;
  
  //std::cerr << "found output " << Quantity << " in PDE " << myPDE_->GetName() << std::endl;

  // create new Coupling Output
  outputTypes_.push_back(outputType);
  outputQuantities_.push_back(quantity);

  CouplingInterface *myInterface = new CouplingInterface;  
  myInterface->region = region;
  myInterface->regionType = regionType;
  myInterface->level = level;

  std::list<Integer> nodesConverted;
  std::vector<Elem*> SD;
  std::list<Integer>::iterator it;
  Integer inode = 0;


   // Get Elements/nodes of coupling region
  switch (regionType)
    {
    case SUBDOMAIN:
      ptGrid_->GetElemSD(SD, region, level);
      ptGrid_->CalcNumberOfNodesInPatch(SD, myInterface->nodes);
      myInterface->size = myInterface->nodes.size();
      myInterface->values.resize(myInterface->nodes.size());
      myInterface->oldValues.resize(myInterface->nodes.size());
 
      
      break;

    case NODES:

      nodesConverted = ptBCs_->GetNodesLevel(region, level);
      myInterface->nodes.resize(nodesConverted.size());
      it = nodesConverted.begin();
      
      inode = 0;
      for (it=nodesConverted.begin(); it != nodesConverted.end(); it++, inode++)
	{
	myInterface->nodes[inode] = *it;
	}
      myInterface->size = myInterface->nodes.size();
      myInterface->values.resize(myInterface->nodes.size());
      myInterface->oldValues.resize(myInterface->nodes.size());

      break;

    case ELEMS1D:
      myInterface->elements = ptBCs_->getEdgesBC(region, level);
      myInterface->size = myInterface->elements.size();
      myInterface->oldValues.resize(myInterface->elements.size());
      
      break;

    case ELEMS2D:
      myInterface->elements = ptBCs_->getFacesBC(region, level);
      myInterface->size = myInterface->elements.size();
      myInterface->values.resize(myInterface->elements.size());
      myInterface->oldValues.resize(myInterface->elements.size());
      
      break;
    }

  outputInterfaces_.push_back(myInterface);
  
  return myInterface;
}




void PDECoupling::SetPDE(BasePDE * aPDE)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::SetPDE" << std::endl;
#endif 
  
  myPDE_ = aPDE;

}

void PDECoupling::SetOutputSize(Integer i, Integer size)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::SetOutputSize" << std::endl;
#endif 
  
  // std::cerr << "entering PDECoupling::SetOutputSize to " << Size << std::endl;
  // std::cerr << "PdeName: " << myPDE_->GetName() << std::endl;

  outputInterfaces_[i]->size = size;
  outputInterfaces_[i]->values.resize(size);
  outputInterfaces_[i]->oldValues.resize(size);

}  

void PDECoupling::SetOutputDim(Integer i, ShortInt dim)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::SetOutputDim" << std::endl;
#endif 
  // std::cerr << "entering PDECoupling::SetOutputDim to " << Dim << std::endl;
  // std::cerr << "PdeName: " << myPDE_->GetName() << std::endl;

  outputInterfaces_[i]->dim = dim;
  outputInterfaces_[i]->values.redim(dim);
  outputInterfaces_[i]->oldValues.redim(dim);
}

std::string PDECoupling::GetPDEName()
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::GetPDEName" << std::endl;
#endif 

  return myPDE_->GetName();
}


Integer PDECoupling::GetNumInputCouplings()
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::GetnumInputCouplings" << std::endl;
#endif 

  return inputInterfaces_.size();
}


Integer PDECoupling::GetNumOutputCouplings()
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::GetnumOutputCouplings" << std::endl;
#endif

  return outputQuantities_.size();
}


CouplingOutputType PDECoupling::Input2OutputType(CouplingInputType inputType)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::Input2OutputType" << std::endl;
#endif 

  switch (inputType)
    {
    case COORD:
      // Coordinate Coupling -> node values needed
      return NODE;
      break;

    case RHS:
      // RHS coupling -> node values needed
      return NODE;
      break;

    case ID_BC:
      // Inhomogeneous Dirichlet BC -> node values needed
      return NODE;
      break;

     case MAT:
      // Material parameter coupling -> element values needed
      return ELEM;
      break;
    }

}



} // end of namespace


