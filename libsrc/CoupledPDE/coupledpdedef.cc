#include "coupledpdedef.hh"
#include "Domain/grid.hh"
#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{


  CoupledPDEDef::CoupledPDEDef(Grid * aptGrid)
  {
    ENTER_FCN( "CoupledPDEDef::CoupledPDEDef" );

    ptGrid_ = aptGrid;

    // Define Ordering of PDEs
    // Here the hardcoded information from coupledPDE.conf is
    // read in and the possible couplings of the different PDEs
    // are defined.
    DefineOrdering();
  }
  

  CoupledPDEDef::~CoupledPDEDef()
  {
    ENTER_FCN( "CoupledPDEDef::~CoupledPDEDef" );

    for (UInt i=0; i<CoupledPDEs_.GetSize(); i++)
      if (CoupledPDEs_[i]) delete CoupledPDEs_[i];
  }
  
  void CoupledPDEDef::CreateCoupling(StdVector<StdPDE*> & OrderedPDEs, 
                                     StdVector<PDECoupling*> & Couplings,
                                     StdVector<StdPDE*> & UnorderedPDEs)
  {
    ENTER_FCN( "CoupledPDEDef::OrderPDEs" );


    bool found = false;
    UInt CoupledPDENumber;
    StdVector<std::string> PDENames;
    OrderedPDEs.Clear();

    //    for (UInt i=0; i<UnorderedPDEs.GetSize(); i++)
    //      std::cerr << "Unorderered PDEs = " << UnorderedPDEs[i]->GetName() << std::endl;

    // iterate over all coupling PDEs to find the 
    // corresponding coupling definition for current set of PDEs
    for (UInt i=0; i<CoupledPDEs_.GetSize(); i++)
      {
        // check if number of PDEs in coupling matches
        if (CoupledPDEs_[i]->GetNumPDEs() == UnorderedPDEs.GetSize())
          {
            CoupledPDEs_[i]->GetNamePDEs(PDENames);
        
            // iterate over all PDEnames in ordered direction
            for (UInt j=0; j<PDENames.GetSize(); j++)
              {
                // iterate over all PDEnames in the vector of unordered PDEs
                for (UInt k=0; k<UnorderedPDEs.GetSize(); k++)
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

    StdVector<std::string> keyVec, attrVec, valVec;
    StdVector<std::string> couplingTerms;


    // iterate over all PDEs specified CoupledPDE
    for (UInt i=0; i<MyCoupledPDE->GetNumPDEs(); i++)
      {
        MyCoupledPDE->GetCouplingType(OrderedPDEs[i]->GetName(), InputType);
        MyCoupledPDE->GetCouplingQuantity(OrderedPDEs[i]->GetName(), InputQuantity);
        MyCoupledPDE->GetCouplingOptionality(OrderedPDEs[i]->GetName(), inputOptionality);
        Couplings[i] = new PDECoupling(ptGrid_);
        Couplings[i]->SetPDE(OrderedPDEs[i]);

        // add all coupling terms of PDE
        for (UInt j=0; j<InputType.GetSize(); j++)

          // if this coupling type is not needed in every coupled simulation
          if (inputOptionality[j])
            {
              keyVec = "couplingList", "iterative", OrderedPDEs[i]->GetName() ,"coupling", "quantity";
              attrVec = "", "", "", "";
              valVec = "", "", "" ,"";
              params->GetList( keyVec, attrVec, valVec, couplingTerms);

              couplingTermsConv.Clear();
              couplingTermsConv.Resize(couplingTerms.GetSize());
              for (UInt k=0; k<couplingTerms.GetSize(); k++)
                String2Enum(couplingTerms[k],couplingTermsConv[k]);

              Boolean found = FALSE;
              for (UInt k=0; k<couplingTermsConv.GetSize(); k++)
                if (couplingTermsConv[k] == InputQuantity[j])
                  found = TRUE;

              if (found) {
                Couplings[i]->RegisterInput(InputType[j], InputQuantity[j]);
              }
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

