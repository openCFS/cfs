#include "pdecoupling.hh"

#include <PDE/basepde.hh>
#include <Domain/grid.hh>
#include <Domain/bcs.hh>
#include <Domain/elem.hh>
#include <list>

namespace CoupledField
{

PDECoupling::PDECoupling(Grid * aptgrid, BCs * aptBCs)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::PDECoupling" << std::endl;
#endif 

  ptGrid_ = aptgrid;
  ptBCs_ = aptBCs;
}

PDECoupling::~PDECoupling()
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::~PDECoupling" << std::endl;
#endif 
  
  for (Integer i=0; i<OutputInterfaces_.size(); i++)
    if (OutputInterfaces_[i]) delete OutputInterfaces_[i];

   for (Integer i=0; i<InputInterfaces_.size(); i++)
    if (InputInterfaces_[i]) delete InputInterfaces_[i];

}

void PDECoupling::RegisterInput(CouplingInputType InType, std::string Quantity)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::RegisterInput" << std::endl;
#endif 
  
  InputTypes_.push_back(InType);
  InputQuantities_.push_back(Quantity);  
  InputInterfaces_.push_back(0);
}

void PDECoupling::AddInput(std::string Quantity, 
			   std::string region, 
			   CouplingRegionType RegionType,
			   Integer level,
			   std::vector<PDECoupling*> & Couplings)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::AddInput" << std::endl;
#endif 

  //std::cerr << "AddInput called with Quantity \"" << Quantity << "\", region \"" << region << "\" regiontype \"" << RegionType << "\" and a vector of Couplings..." << std::endl;

  CouplingInterface *MyInterface = 0; 
  Integer MyNum = -1;

  // search matching  Quantity
  for (Integer i=0; i<InputQuantities_.size(); i++)
    {
      if (InputQuantities_[i] == Quantity)
	MyNum = i; 
    }

   if (MyNum == -1)
    {
      std::string ErrMsg = "Qantity \'" + Quantity + "\' not registered for PDE \'" + myPDE_->GetName()+ "\'";
      Error(ErrMsg.c_str(),__FILE__,__LINE__);
    }

 
  // Add coupling input as output to pde, which can calculate specified quantity
  CouplingOutputType MyOutputType = Input2OutputType(InputTypes_[MyNum]);

  Integer i=0;
  while (MyInterface == 0 && i<Couplings.size())
    {
      MyInterface = Couplings[i++]->AddOutput(MyOutputType, Quantity, region, RegionType, level);
    }
  
  // If no pde has the specified quantity as output
  if (MyInterface == 0)
    {
      std::string ErrMsg = "Qantity \'" + Quantity + "\' can not be calculated with current set of PDEs";
      Error(ErrMsg.c_str(),__FILE__,__LINE__);
    }
  
  InputInterfaces_[MyNum] = MyInterface;

}

PDECoupling::CouplingInterface* PDECoupling::AddOutput(CouplingOutputType OutputType, 
						       std::string Quantity, 
						       std::string region,
						       CouplingRegionType RegionType,
						       Integer level)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::AddOutput" << std::endl;
#endif 

  // search if output exists already
  for (Integer i=0; i<OutputTypes_.size(); i++)
    if (OutputTypes_[i] == OutputType 
	&& OutputQuantities_[i] == Quantity
	&& OutputInterfaces_[i]->Region == region)
      return OutputInterfaces_[i];
  
  
  // if not find out if PDE can calculate desired output
  if (!myPDE_->HasOutput(Quantity))
    return 0;
  
  // create new Coupling Output
  OutputTypes_.push_back(OutputType);
  OutputQuantities_.push_back(Quantity);
  CouplingInterface *MyInterface = new CouplingInterface;  
  MyInterface->Region = region;
  MyInterface->RegionType = RegionType;
  MyInterface->level = level;
  std::list<Integer> nodesConverted;
  std::vector<Elem*> SD;
  std::list<Integer>::iterator it;
  Integer inode = 0;


   // Get Elements/nodes of coupling region
  switch (RegionType)
    {
    case SUBDOMAIN:
      ptGrid_->GetElemSD(SD, region, level);
      ptGrid_->CalcNumberOfNodesInPatch(SD, MyInterface->Nodes);
      MyInterface->Values.resize(MyInterface->Nodes.size());
      MyInterface->Size = MyInterface->Nodes.size();
      
      break;

    case NODES:

      nodesConverted = ptBCs_->GetNodesLevel(region, level);
      MyInterface->Nodes.resize(nodesConverted.size());
      it = nodesConverted.begin();
      
      inode = 0;
      for (it=nodesConverted.begin(); it != nodesConverted.end(); it++, inode++)
	{
	MyInterface->Nodes[inode] = *it;
	}

      MyInterface->Values.resize(MyInterface->Nodes.size());
      MyInterface->Size = MyInterface->Nodes.size();
      
      break;

    case ELEMS1D:
      MyInterface->Elements = ptBCs_->getEdgesBC(region, level);

      MyInterface->Values.resize(MyInterface->Elements.size());
      MyInterface->Size = MyInterface->Elements.size();

      break;

    case ELEMS2D:
      MyInterface->Elements = ptBCs_->getFacesBC(region, level);

      MyInterface->Values.resize(MyInterface->Elements.size());
      MyInterface->Size = MyInterface->Elements.size();
      
      break;
    }

  OutputInterfaces_.push_back(MyInterface);
  
  return MyInterface;
}




void PDECoupling::SetPDE(BasePDE * aPDE)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::SetPDE" << std::endl;
#endif 
  
  myPDE_ = aPDE;

}

void PDECoupling::SetOutputSize(Integer i, Integer Size)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::SetOutputSize" << std::endl;
#endif 
  
  std::cerr << "entering PDECoupling::SetOutputSize to " << Size << std::endl;
  std::cerr << "PdeName: " << myPDE_->GetName() << std::endl;
  OutputInterfaces_[i]->Size = Size;
  OutputInterfaces_[i]->Values.resize(Size);

}  

void PDECoupling::SetOutputDim(Integer i, ShortInt Dim)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::SetOutputDim" << std::endl;
#endif 
  std::cerr << "entering PDECoupling::SetOutputDim to " << Dim << std::endl;
  std::cerr << "PdeName: " << myPDE_->GetName() << std::endl;
  OutputInterfaces_[i]->Dim = Dim;
  OutputInterfaces_[i]->Values.redim(Dim);

}

Integer PDECoupling::GetNumInputCouplings()
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::GetnumInputCouplings" << std::endl;
#endif 

  return InputInterfaces_.size();
}


Integer PDECoupling::GetNumOutputCouplings()
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::GetnumOutputCouplings" << std::endl;
#endif

  return OutputInterfaces_.size();
}


Double PDECoupling::CalcNorm(std::string type)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::CalcNorm" << std::endl;
#endif 

}


CouplingOutputType PDECoupling::Input2OutputType(CouplingInputType InputType)
{
#ifdef TRACE
  (*trace) << "entering PDECoupling::Input2OutputType" << std::endl;
#endif 

  switch (InputType)
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


