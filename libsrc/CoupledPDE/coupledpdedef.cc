#include "coupledpdedef.hh"
#include <Domain/grid.hh>
#include <Domain/bcs.hh>

namespace CoupledField
{


CoupledPDEDef::CoupledPDEDef(Grid * aptGrid, BCs * aptBCs)
{
#ifdef TRACE
  (*trace) << "entering  CoupledPDEDef::CoupledPDEDef" << std::endl;
#endif 

  ptGrid_ = aptGrid;
  ptBCs_ = aptBCs;

  // Define Ordering of PDEs
  DefineOrdering();
}
  

CoupledPDEDef::~CoupledPDEDef()
{
#ifdef TRACE
  (*trace) << "entering  CoupledPDEDef::~CoupledPDEDef" << std::endl;
#endif 

  for (Integer i=0; i<CoupledPDEs_.size(); i++)
    if (CoupledPDEs_[i]) delete CoupledPDEs_[i];
}
  
void CoupledPDEDef::CreateCoupling(std::vector<BasePDE*> & OrderedPDEs, 
				   std::vector<PDECoupling*> & Couplings,
				   std::vector<BasePDE*> & UnorderedPDEs)
{
#ifdef TRACE
  (*trace) << "entering  CoupledPDEDef::OrderPDEs" << std::endl;
#endif   

   bool found = false;
   Integer CoupledPDENumber;
   std::vector<std::string> PDENames;
   OrderedPDEs.clear();

   // iterate over all coupling PDEs
   //std::cerr << "CoupledPDEs_.size() = " << CoupledPDEs_.size() << std::endl;
   //std::cerr << "UnorderedPDEs.size() = " << UnorderedPDEs.size() << std::endl;

   for (Integer i=0; i<CoupledPDEs_.size(); i++)
     {
       
      // check if number of PDEs in coupling matches
      if (CoupledPDEs_[i]->GetNumPDEs() == UnorderedPDEs.size())
	{
	  //std::cerr << "CoupledPDEs[" << i << ".NumPDE == UnorderedPDEs.size()" << std::endl;
	  CoupledPDEs_[i]->GetNamePDEs(PDENames);
	
	  //std::cerr << "PDENames.size() = " << PDENames.size() << std::endl;
	  // iterate over all PDEnames in ordered direction
	  for (Integer j=0; j<PDENames.size(); j++)
	    {
	      //std::cerr << "PDENames[" << j << "] = " << PDENames[j] << std::endl;
	      // iterate over all PDEnames in the vector of unordered PDEs
	      for (Integer k=0; k<UnorderedPDEs.size(); k++)
		if (PDENames[j] == UnorderedPDEs[k]->GetName())
		  OrderedPDEs.push_back(UnorderedPDEs[k]);
	    }
	}
      
      // check if all PDEs could be assigned
      if (OrderedPDEs.size() == CoupledPDEs_[i]->GetNumPDEs())
	{
	  found = true;
	  CoupledPDENumber = i;
	  break;
	} 
      else
	{
	  OrderedPDEs.clear();
	}
      
    }
  
  if (found != true)
    Error("Coupling for current set of PDEs ist not defined!",__FILE__,__LINE__);

  
  // Create Coupling objects
  Couplings.clear();
  Definition * MyCoupledPDE = CoupledPDEs_[CoupledPDENumber];					   
  std::vector<CouplingInputType>  InputType;
  std::vector<std::string> InputQuantity;
  Couplings.resize(MyCoupledPDE->GetNumPDEs());

  // iterate over all PDEs specified CoupledPDE
  for (Integer i=0; i<MyCoupledPDE->GetNumPDEs(); i++)
    {
      MyCoupledPDE->GetCouplingType(OrderedPDEs[i]->GetName(), InputType);
      MyCoupledPDE->GetCouplingQuantity(OrderedPDEs[i]->GetName(), InputQuantity);
      Couplings[i] = new PDECoupling(ptGrid_, ptBCs_);
      //std::cerr << "OrderedPDEs[i]->Name() = " << OrderedPDEs[i]->GetName() << std::endl;
      Couplings[i]->SetPDE(OrderedPDEs[i]);

      // add all coupling terms of PDE
      for (Integer j=0; j<InputType.size(); j++)
	{
	  //std::cerr << "Adding to " << OrderedPDEs[i]->GetName() << " " << InputQuantity[j] << std::endl;
	    Couplings[i]->RegisterInput(InputType[j], InputQuantity[j]);  
	}
    }
  
  //std::cerr << "Finished Create Coupling" << std::endl;

}

void CoupledPDEDef::DefineOrdering()
{
#ifdef TRACE
  (*trace) << "entering  CoupledPDEDef::DefineOrdering" << std::endl;
#endif 


#include <CoupledPDE/coupledPDE.conf>
}

Definition::Definition()
{
#ifdef TRACE
  (*trace) << "entering  Definition::Definition" << std::endl;
#endif 
}

Definition::~Definition()
{
#ifdef TRACE
  (*trace) << "entering  Definition::~Definition" << std::endl;
#endif
}

void Definition::AddPDE(std::string PDEName)
{
#ifdef TRACE
  (*trace) << "entering  Definition::AddPDE" << std::endl;
#endif 

  PDEs_.push_back(PDEName);
  NumPDEs_ = PDEs_.size();
}

void Definition::AddInputCoupling(std::string PDEName, CouplingInputType InType, std::string Quantity)
{
#ifdef TRACE
  (*trace) << "entering  Definition::AddCoupling" << std::endl;
#endif 

  InputCouplingTypes_[PDEName].push_back(InType);
  InputCouplingQuantities_[PDEName].push_back(Quantity);
}


} // end of namespace

