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
  // Here the hardcoded information from coupledPDE.conf is
  // read in and the possible couplings of the different PDEs
  // are defined.
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

   // iterate over all coupling PDEs to find the 
   // corresponding coupling definition for current set of PDEs
   for (Integer i=0; i<CoupledPDEs_.size(); i++)
     {
       
      // check if number of PDEs in coupling matches
      if (CoupledPDEs_[i]->GetNumPDEs() == UnorderedPDEs.size())
	{
	  CoupledPDEs_[i]->GetNamePDEs(PDENames);
	
	  // iterate over all PDEnames in ordered direction
	  for (Integer j=0; j<PDENames.size(); j++)
	    {
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
  std::vector<Boolean> inputOptionality;
  Couplings.resize(MyCoupledPDE->GetNumPDEs());

  std::vector<std::string> couplingTerms;

  // iterate over all PDEs specified CoupledPDE
  for (Integer i=0; i<MyCoupledPDE->GetNumPDEs(); i++)
    {
      MyCoupledPDE->GetCouplingType(OrderedPDEs[i]->GetName(), InputType);
      MyCoupledPDE->GetCouplingQuantity(OrderedPDEs[i]->GetName(), InputQuantity);
      MyCoupledPDE->GetCouplingOptionality(OrderedPDEs[i]->GetName(), inputOptionality);
      Couplings[i] = new PDECoupling(ptGrid_, ptBCs_);
      Couplings[i]->SetPDE(OrderedPDEs[i]);

      // add all coupling terms of PDE
      for (Integer j=0; j<InputType.size(); j++)

	// if this coupling type is not needed every coupled simulation
	if (inputOptionality[j])
	  {
	    conf->getliststr("input_coupling_terms", couplingTerms, OrderedPDEs[i]->GetName());
	    
	    Boolean found = FALSE;
	    
	    for (Integer k=0; k<couplingTerms.size(); k++)
	      if (couplingTerms[k] == InputQuantity[j])
		found = TRUE;

	    if (found)
	      Couplings[i]->RegisterInput(InputType[j], InputQuantity[j]);
	  }
	else
	  Couplings[i]->RegisterInput(InputType[j], InputQuantity[j]);
    }
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

void Definition::AddInputCoupling(std::string PDEName, 
				  CouplingInputType InType, 
				  std::string Quantity,
				  Boolean optionalCoupling) //"optionalCoupling" is by default FALSE
{
#ifdef TRACE
  (*trace) << "entering  Definition::AddCoupling" << std::endl;
#endif 

  InputCouplingTypes_[PDEName].push_back(InType);
  InputCouplingQuantities_[PDEName].push_back(Quantity);
  optionalCoupling_[PDEName].push_back(optionalCoupling);  
}


} // end of namespace

