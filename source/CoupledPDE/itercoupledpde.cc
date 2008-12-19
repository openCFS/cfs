// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "itercoupledpde.hh"

#include "pdecoupling.hh"

#include "DataInOut/WriteInfo.hh"
#include "PDE/SinglePDE.hh"
#include "Driver/iterSolveStep.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"

namespace CoupledField
{

  IterCoupledPDE::IterCoupledPDE(StdVector<StdPDE*> & PDEs,
                                 StdVector<SinglePDE*> & singlePDEs,
                                 StdVector<PDECoupling*> & Couplings,
                                 ParamNode* paramNode )
    : BasePDE( paramNode )
  {

    PDEs_       = PDEs;
    singlePDEs_ = singlePDEs;
    Couplings_  = Couplings;
    myParam_ = paramNode;

    NumPDEs_ = PDEs.GetSize();
    solveStep_ = NULL;

    // Concatenate PDE name strings for output in info-file
    pdename_ = "CoupledPDE: ";
    for (UInt actPDE=0; actPDE < PDEs.GetSize()-1; actPDE++)
      pdename_ += PDEs[actPDE] -> GetName() + "+";
    pdename_ += PDEs[PDEs.GetSize()-1] -> GetName();

    // fetch "nonLinear" node
    ParamNode * nonLinNode = myParam_->Get("nonLinear", false );

    // get maximum number of iterations (optional)
    maxiter_ = 100;
    if( nonLinNode )
      nonLinNode->Get( "maxNumIters", maxiter_, false );

    // query logging flag
    nonLinLogging_ = true;
    if( nonLinNode )
      nonLinNode->Get( "logging", nonLinLogging_, false );

  }


  IterCoupledPDE::~IterCoupledPDE()
  {

    // delete solveStep-object
    delete solveStep_;

    // delete coupling objects
    for ( UInt i = 0; i < Couplings_.GetSize(); i++ ) {
      delete Couplings_[i];
    }
    Couplings_.Clear();
  }


  void IterCoupledPDE::InitCoupling()
  {

    StdVector<std::string> quantities;
    StdVector<std::string> interfaceTypes;
    StdVector<std::string> interfaceNames;
    StdVector<std::string> stopCritQuantities;
    StdVector<std::string> neighbourRegions;
    std::string normtype;
    std::string errMsg;
    Double epsilon;

    // vectors containing each quantity only one
    // time
    StdVector<std::string> quantitiesSorted;
    StdVector<std::string> interfaceTypesSorted;
    StdVector<StdVector<std::string> > interfaceNamesSorted;
    StdVector<std::string> nameVectorAux;


    // Iterate over all PDEs
    for ( UInt iPDE = 0; iPDE < singlePDEs_.GetSize(); iPDE++ ) {

      quantities.Clear();
      interfaceTypes.Clear();
      interfaceNames.Clear();
      quantitiesSorted.Clear();
      interfaceTypesSorted.Clear();
      interfaceNamesSorted.Clear();

      // get parameter nodes of current pde and look for the one
      // with the same name as the curent searched-for pde
      ParamNode * pdeNode = NULL;
      for( UInt i = 0; i < myParam_->GetChildren().GetSize(); i++ ) {
        // check, id node is "nonLinear"
        if( (myParam_->GetChildren())[i]->GetName() == "nonLinear")
          continue;

        // fetch names of related pdes
        StdVector<ParamNode*> pdeNodes = (myParam_->GetChildren()[i])->GetChildren();

        for(UInt iNode = 0; iNode < pdeNodes.GetSize(); iNode++ ) {
          std::string name = pdeNodes[iNode]->GetName();
          if( name == singlePDEs_[iPDE]->GetName() ) {
            pdeNode = pdeNodes[iNode];
            break;
          }
        }

      }
      if (!pdeNode)
        EXCEPTION( "PDE '" << singlePDEs_[iPDE]->GetName() << "' could "
                   << " not be found in the iterative coupling section");


      // get coupling entities
      StdVector<ParamNode*> couplNodes = pdeNode->GetList("coupling");

      quantities.Clear();
      interfaceTypes.Clear();
      interfaceNames.Clear();

      // iterate over all children nodes
      for( UInt i = 0; i < couplNodes.GetSize(); i++ ) {

        // get quantity
        quantities.Push_back( couplNodes[i]->Get("quantity")->AsString() );

        // get type
        interfaceTypes.Push_back( couplNodes[i]->Get("type")->AsString() );

        // get name
        interfaceNames.Push_back( couplNodes[i]->Get("name")->AsString() );
      }

      std::string typeAux;
      std::string quantityTemp;
      std::string nameAux;

      // Check if for one quantity different kinds of
      // interfaces were specified
      for (UInt i=0; i<quantities.GetSize(); i++) {

        typeAux = interfaceTypes[i];
        quantityTemp = quantities[i];

        for (UInt j=0; j<quantities.GetSize(); j++) {
          if (quantityTemp == quantities[j]) {
            if (typeAux != interfaceTypes[j]) {
              EXCEPTION( "Per coupling quantity only one type (node, region"
                         << ", surface) of interface is allowed. \n You specified "
                         << "for quantity '" << quantityTemp << "' as interface '"
                         << typeAux << "' and '"  << interfaceTypes[j] << "."
                         << "Please correct the parameter file!"; );
            }
          }
        }
      }

      // Now merge for each coupling quantity all interface names
      // of the same type (node, interface, region)
      bool found = false;
      for (UInt iQuant=0; iQuant<quantities.GetSize(); iQuant++) {
        typeAux = interfaceTypes[iQuant];
        quantityTemp = quantities[iQuant];
        nameAux = interfaceNames[iQuant];
        found = false;
        for (UInt iSorted=0; iSorted<quantitiesSorted.GetSize(); iSorted++){
          if(quantityTemp == quantitiesSorted[iSorted]) {
            interfaceNamesSorted[iSorted].Push_back(nameAux);
            found = true;
          }
        }
        if (!found) {
          quantitiesSorted.Push_back(quantityTemp);
          interfaceTypesSorted.Push_back(typeAux);
          nameVectorAux.Clear();
          nameVectorAux.Push_back(nameAux);
          interfaceNamesSorted.Push_back(nameVectorAux);
        }
      }

      NormType normTypeAux;
      CouplingRegionType regionTypeAux;
      SolutionType quantityAux;

      // Get list for all quantities, which are listed
      // in section "nonLinear" and have a stopping
      // criterion
      stopCritQuantities.Clear();
      ParamNode * nonLinNode = myParam_->Get("nonLinear", false);
      StdVector<ParamNode*> stopNodes;
      if( nonLinNode ) {
        stopNodes = nonLinNode->GetList("stopCrit");
        stopCritQuantities.Clear();
        for( UInt i = 0; i < stopNodes.GetSize(); i++ )
          stopCritQuantities.Push_back( stopNodes[i]->Get("quantity")->AsString() );
      }

      // Get for each quantity the according stopping Criteria and normtype
      for( UInt iQuant = 0; iQuant < quantitiesSorted.GetSize(); iQuant++ ){

        epsilon = 0.0;
        normtype.clear();

        String2Enum(quantitiesSorted[iQuant], quantityAux);

        Integer index = stopCritQuantities.Find(quantitiesSorted[iQuant]);

        // check, if a stopping criterion was defined for current quantity
        if ( index != -1 ) {
          stopNodes[index]->Get( "value", epsilon ); // stopping criterion
          stopNodes[index]->Get( "l2Norm", normtype ); // type of norm used

          // if quantity is elecForceVWP or magForceVWP, get neighbouring region
          neighbourRegions.Clear();
          if (quantityAux == MAG_FORCE_VWP ||
              quantityAux == ELEC_FORCE_VWP) {
            StdVector<ParamNode*> regionNodes =
              pdeNode->GetList("coupling", "quantity", quantitiesSorted[iQuant]);

            for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
              neighbourRegions.Push_back( regionNodes[i]->Get("neighbourRegion")->AsString() );
            }
          }

        }
        else {
          epsilon = 0.0;
          normtype = "no";
        }

        String2Enum(normtype, normTypeAux);
        String2Enum(interfaceTypesSorted[iQuant], regionTypeAux);



        // register the interface at the according coupling-object
        Couplings_[iPDE]->AddInput(quantityAux,
                                   interfaceNamesSorted[iQuant],
                                   regionTypeAux, neighbourRegions,
                                   epsilon, normTypeAux, Couplings_);
        norms_.Push_back(1.0);
      }
    }

    // Initialize each PDEs coupling terms
    for ( UInt i = 0; i < Couplings_.GetSize(); i++ ) {
      Couplings_[i]->GetPDE()->InitCoupling(Couplings_[i]);
    }

    ToInfo(info->Get(InfoNode::HEADER)->Get("coupling"));

    // create solve step object
    solveStep_ = new IterSolveStep( *this );
  }


  void IterCoupledPDE::WriteRestart()
  {

    for (UInt actPDE=0; actPDE < PDEs_.GetSize(); actPDE++)
      PDEs_[actPDE]->WriteRestart( );
  }


  void IterCoupledPDE::ReadRestart(UInt &startStep)
  {

    for (UInt actPDE=0; actPDE < PDEs_.GetSize(); actPDE++)
      PDEs_[actPDE]->ReadRestart(startStep);
  }


  void IterCoupledPDE::WriteResultsInFile(const UInt kstep,
                                          const Double asteptime )
  {

    for (UInt i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->WriteResultsInFile(kstep, asteptime );
  }


  void IterCoupledPDE::ToInfo(InfoNode* base)
  {
    CFSVector *val;
    StdVector<UInt> * nodes;
    StdVector<std::string> couplingRegions;
    std::string couplingTypeAux, quantityAux, regionTypeAux, normTypeAux;

    for (UInt ipde=0; ipde<PDEs_.GetSize(); ipde++)
    {
      InfoNode* pde = base->Get(Couplings_[ipde]->GetPDE()->GetName());

      // Show InputCouplings
      InfoNode* ic = pde->Get("inputCouplings");
      for (UInt i=0; i<Couplings_[ipde]->GetNumInputCouplings(); i++)
      {
        InfoNode* c = ic->Get("coupling", InfoNode::APPEND);
        Couplings_[ipde]->GetInputNodes(i, nodes);
        Couplings_[ipde]->GetInputValues(i,val);

        // Convert enum-types in strings
        Enum2String(Couplings_[ipde]->GetInputType(i), couplingTypeAux);
        Enum2String(Couplings_[ipde]->GetInputQuantity(i), quantityAux);
        Enum2String(Couplings_[ipde]->GetInputRegionType(i), regionTypeAux);
        Enum2String(Couplings_[ipde]->GetInputNormType(i), normTypeAux);

        c->Get("inputCoupling")->SetValue(i+1);
        c->Get("couplingType")->SetValue(couplingTypeAux);
        c->Get("inputQuantity")->SetValue(quantityAux);

        InfoNode* list = c->Get("regions");
        Couplings_[ipde]->GetInputRegions(i, couplingRegions);
        for (UInt j=0; j<couplingRegions.GetSize()-1; j++)
          list->Get("region", InfoNode::APPEND)->Get("name")->SetValue(couplingRegions[j]);

        c->Get("regionType")->SetValue(regionTypeAux);
        c->Get("numberInputCouplingValues")->SetValue(val->GetSize());
        c->Get("inputValueDof")->SetValue(Couplings_[ipde]->GetInputDof(i));
        c->Get("numberInputNodes")->SetValue(Couplings_[ipde]->GetInputNumNodes(i));
        c->Get("numberInputElems")->SetValue(Couplings_[ipde]->GetInputNumElems(i));
        c->Get("normType")->SetValue(normTypeAux);
        c->Get("tolerance")->SetValue(Couplings_[ipde]->GetInputEpsilon(i));
      }

      couplingRegions.Clear();
      // Show OutputCouplings
      nodes = 0;

      InfoNode* oc = pde->Get("outputCouplings");
      for (UInt i=0; i<Couplings_[ipde]->GetNumOutputCouplings(); i++)
      {
        InfoNode* c = oc->Get("coupling", InfoNode::APPEND);
        Couplings_[ipde]->GetOutputNodes(i, nodes);
        Couplings_[ipde]->GetOutputValues(i,val);

        // Convert enum-types in strings
        Enum2String(Couplings_[ipde]->GetOutputType(i), couplingTypeAux);
        Enum2String(Couplings_[ipde]->GetOutputQuantity(i), quantityAux);
        Enum2String(Couplings_[ipde]->GetOutputRegionType(i), regionTypeAux);
        Enum2String(Couplings_[ipde]->GetOutputNormType(i), normTypeAux);

        c->Get("outputCoupling")->SetValue(i+1);
        c->Get("couplingType")->SetValue(couplingTypeAux);
        c->Get("inputQuantity")->SetValue(quantityAux);

        InfoNode* list = c->Get("regions");
        Couplings_[ipde]->GetOutputRegions(i, couplingRegions);
        for (UInt j=0; j<couplingRegions.GetSize()-1; j++)
          list->Get("region", InfoNode::APPEND)->Get("name")->SetValue(couplingRegions[j]);

        c->Get("regionType")->SetValue(regionTypeAux);
        c->Get("numberOutputCouplingValues")->SetValue(val->GetSize());
        c->Get("outputValueDof")->SetValue(Couplings_[ipde]->GetOutputDof(i));
        c->Get("numberOutputNodes")->SetValue(Couplings_[ipde]->GetOutputNumNodes(i));
        c->Get("numberOutputElems")->SetValue(Couplings_[ipde]->GetOutputNumElems(i));
        c->Get("normType")->SetValue(normTypeAux);
        c->Get("tolerance")->SetValue(Couplings_[ipde]->GetOutputNumElems(i));
      }
    }
  }


  void IterCoupledPDE::WriteGeneralPDEdefines()
  {

    for (UInt i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->WriteGeneralPDEdefines();
  }


  BaseSolveStep * IterCoupledPDE::GetSolveStep()
  {
    return solveStep_;
  }


} // end of namespace
