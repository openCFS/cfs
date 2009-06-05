// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "coupledpdedef.hh"
#include "Domain/grid.hh"
#include "General/environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{


  CoupledPDEDef::CoupledPDEDef(Grid * aptGrid)
  {

    ptGrid_ = aptGrid;

    // Define Ordering of PDEs
    // Here the hardcoded information from coupledPDE.conf is
    // read in and the possible couplings of the different PDEs
    // are defined.
    DefineOrdering();
  }
  

  CoupledPDEDef::~CoupledPDEDef()
  {

    for (UInt i=0; i<CoupledPDEs_.GetSize(); i++)
      if (CoupledPDEs_[i]) delete CoupledPDEs_[i];
  }
  
  void CoupledPDEDef::CreateCoupling(StdVector<StdPDE*> & OrderedPDEs, 
                                     StdVector<PDECoupling*> & Couplings,
                                     StdVector<StdPDE*> & UnorderedPDEs,
                                     ParamNode * iterCoupledNode )
  {


    bool found = false;
    UInt CoupledPDENumber = 0;
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
      EXCEPTION("Coupling for current set of PDEs ist not defined!");

    // Create Coupling objects
    Couplings.Clear();
    Definition * MyCoupledPDE = CoupledPDEs_[CoupledPDENumber];                                      
    StdVector<CouplingInputType>  InputType;
    StdVector<SolutionType> InputQuantity;
    StdVector<SolutionType> couplingTermsConv;
    StdVector<bool> inputOptionality;
    Couplings.Resize(MyCoupledPDE->GetNumPDEs());
    Couplings.Init( NULL );

    StdVector<std::string> couplingTerms;


    // iterate over all PDEs specified CoupledPDE
    for (UInt i=0; i<MyCoupledPDE->GetNumPDEs(); i++) {
      MyCoupledPDE->GetCouplingType(OrderedPDEs[i]->GetName(), InputType);
      MyCoupledPDE->GetCouplingQuantity(OrderedPDEs[i]->GetName(), InputQuantity);
      MyCoupledPDE->GetCouplingOptionality(OrderedPDEs[i]->GetName(), inputOptionality);
      Couplings[i] = new PDECoupling(ptGrid_);
      Couplings[i]->SetPDE(OrderedPDEs[i]);
      
      // add all coupling terms of PDE
      for (UInt j=0; j<InputType.GetSize(); j++) {
        
        // if this coupling type is not needed in every coupled simulation
        if (inputOptionality[j]) {
   
          // look for correct pde
          std::string pdeName = OrderedPDEs[i]->GetName();
          ParamNode * pdeNode = NULL;
          for( UInt iNode = 0; iNode < iterCoupledNode->GetChildren().GetSize(); iNode++ ) {
            // check, id node is "nonLinear"
            if( (iterCoupledNode->GetChildren())[iNode]->GetName() == "nonLinear") 
              continue;
            
            if( (iterCoupledNode->GetChildren())[iNode]->GetChildren()[1]->GetName() 
                == pdeName ) {
              pdeNode = (iterCoupledNode->GetChildren())[iNode]->GetChildren()[1];
              break;
            } else if( iterCoupledNode->GetChildren()[iNode]->GetChildren()[2]->GetName() 
                       == pdeName )  {
              pdeNode = (iterCoupledNode->GetChildren())[iNode]->GetChildren()[2];
              break;
            }
          }
          if (!pdeNode) 
            break;

          // create vector with paramnode for each coupling region
          StdVector<ParamNode*> couplNodes = pdeNode->GetList( "coupling");

          couplingTermsConv.Clear();
          couplingTermsConv.Resize(couplNodes.GetSize());
          for (UInt k = 0; k < couplNodes.GetSize(); k++) {
            couplingTermsConv[k] = SolutionTypeEnum.Parse(couplNodes[k]->Get("quantity"));
          }
          
          bool found = false;
          for (UInt k=0; k<couplingTermsConv.GetSize(); k++)
            if (couplingTermsConv[k] == InputQuantity[j])
              found = true;
          
          if (found) {
            Couplings[i]->RegisterInput(InputType[j], InputQuantity[j]);
          }
        }
        else
          Couplings[i]->RegisterInput(InputType[j], InputQuantity[j]);
      }
    }
  }
  



  void CoupledPDEDef::DefineOrdering()
  {

#include <CoupledPDE/coupledPDE.conf>
  }

  Definition::Definition()
  {
  }

  Definition::~Definition()
  {
  }

  void Definition::AddPDE(std::string PDEName)
  {

    PDEs_.Push_back(PDEName);
    NumPDEs_ = PDEs_.GetSize();
  }

  void Definition::AddInputCoupling(std::string PDEName, 
                                    CouplingInputType InType, 
                                    SolutionType Quantity,
                                    bool optionalCoupling) //"optionalCoupling" is by default false
  {

    InputCouplingTypes_[PDEName].Push_back(InType);
    InputCouplingQuantities_[PDEName].Push_back(Quantity);
    optionalCoupling_[PDEName].Push_back(optionalCoupling);  
  }


} // end of namespace

