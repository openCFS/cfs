#include "coupledpdedef.hh"
#include "Domain/grid.hh"
#include "Domain/bcs.hh"
#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{


CoupledPDEDef::CoupledPDEDef(Grid * aptGrid, BCs * aptBCs)
{
  ENTER_FCN( "CoupledPDEDef::CoupledPDEDef" );

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
  ENTER_FCN( "CoupledPDEDef::~CoupledPDEDef" );

  for (Integer i=0; i<CoupledPDEs_.GetSize(); i++)
    if (CoupledPDEs_[i]) delete CoupledPDEs_[i];
}
  
void CoupledPDEDef::CreateCoupling(StdVector<BasePDE*> & OrderedPDEs, 
				   StdVector<PDECoupling*> & Couplings,
				   StdVector<BasePDE*> & UnorderedPDEs)
{
  ENTER_FCN( "CoupledPDEDef::OrderPDEs" );


   bool found = false;
   Integer CoupledPDENumber;
   StdVector<std::string> PDENames;
   OrderedPDEs.Clear();

//    for (Integer i=0; i<UnorderedPDEs.GetSize(); i++)
//      std::cerr << "Unorderered PDEs = " << UnorderedPDEs[i]->GetName() << std::endl;

   // iterate over all coupling PDEs to find the 
   // corresponding coupling definition for current set of PDEs
   for (Integer i=0; i<CoupledPDEs_.GetSize(); i++)
     {
      // check if number of PDEs in coupling matches
      if (CoupledPDEs_[i]->GetNumPDEs() == UnorderedPDEs.GetSize())
	{
	  CoupledPDEs_[i]->GetNamePDEs(PDENames);
	
	  // iterate over all PDEnames in ordered direction
	  for (Integer j=0; j<PDENames.GetSize(); j++)
	    {
	      // iterate over all PDEnames in the vector of unordered PDEs
	      for (Integer k=0; k<UnorderedPDEs.GetSize(); k++)
		if (PDENames[j] == UnorderedPDEs[k]->GetName())
		  OrderedPDEs.Push_back(UnorderedPDEs[k]);
	    }
	}
      
      // check if all PDEs could be assigned
      if (OrderedPDEs.GetSize() == CoupledPDEs_[i]->GetNumPDEs())
	{
	  found = true;
	  CoupledPDENumber = i;
	  break;
	} 
      else
	{
	  OrderedPDEs.Clear();
	}
      
    }
  
  if (found != true)
    Error("Coupling for current set of PDEs ist not defined!",__FILE__,__LINE__);

  
  // Create Coupling objects
  Couplings.Clear();
  Definition * MyCoupledPDE = CoupledPDEs_[CoupledPDENumber];					   
  StdVector<CouplingInputType>  InputType;
  StdVector<SolutionType> InputQuantity;
  StdVector<SolutionType> couplingTermsConv;
  StdVector<Boolean> inputOptionality;
  Couplings.Resize(MyCoupledPDE->GetNumPDEs());

  StdVector<std::string> couplingTerms;


  // iterate over all PDEs specified CoupledPDE
  for (Integer i=0; i<MyCoupledPDE->GetNumPDEs(); i++)
    {
      MyCoupledPDE->GetCouplingType(OrderedPDEs[i]->GetName(), InputType);
      MyCoupledPDE->GetCouplingQuantity(OrderedPDEs[i]->GetName(), InputQuantity);
      MyCoupledPDE->GetCouplingOptionality(OrderedPDEs[i]->GetName(), inputOptionality);
      Couplings[i] = new PDECoupling(ptGrid_, ptBCs_);
      Couplings[i]->SetPDE(OrderedPDEs[i]);

      // add all coupling terms of PDE
      for (Integer j=0; j<InputType.GetSize(); j++)

	// if this coupling type is not needed in every coupled simulation
	if (inputOptionality[j])
	  {
	    params->GetList( "quantity", couplingTerms, "couplingList", "coupling");

	    couplingTermsConv.Clear();
	    couplingTermsConv.Resize(couplingTerms.GetSize());
	    for (Integer k=0; k<couplingTerms.GetSize(); k++)
		 String2Enum(couplingTerms[k],couplingTermsConv[i]);

	    Boolean found = FALSE;
	    
	    for (Integer k=0; k<couplingTermsConv.GetSize(); k++)
	      if (couplingTermsConv[k] == InputQuantity[j])
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
  ENTER_FCN( "CoupledPDEDef::DefineOrdering" );

#include <CoupledPDE/coupledPDE.conf>
}

Definition::Definition()
{
  ENTER_FCN ( "Definition::Definition" );
}

Definition::~Definition()
{
  ENTER_FCN( "Definition::~Definition" );
}

void Definition::AddPDE(std::string PDEName)
{
  ENTER_FCN( "Definition::AddPDE" );

  PDEs_.Push_back(PDEName);
  NumPDEs_ = PDEs_.GetSize();
}

void Definition::AddInputCoupling(std::string PDEName, 
				  CouplingInputType InType, 
				  SolutionType Quantity,
				  Boolean optionalCoupling) //"optionalCoupling" is by default FALSE
{
  ENTER_FCN( "Definition::AddCoupling" );

  InputCouplingTypes_[PDEName].Push_back(InType);
  InputCouplingQuantities_[PDEName].Push_back(Quantity);
  optionalCoupling_[PDEName].Push_back(optionalCoupling);  
}


} // end of namespace

