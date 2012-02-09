// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <iostream>
#include <map>
#include <utility>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "Domain/surfElem.hh"
#include "General/Enum.hh"
#include "General/exception.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/vector.hh"
#include "PDE/StdPDE.hh"
#include "couplingmemento.hh"
#include "pdecoupling.hh"

// header for logging
#include "DataInOut/Logging/cfslog.hh"



namespace CoupledField
{
// declare logging stream
DECLARE_LOG(pdecpl)
DEFINE_LOG(pdecpl, "pdecpl")

  PDECoupling::CouplingInterface::CouplingInterface()
  {
  
    dof = 0;
    numNodes = 0;    
    numElems = 0;
    epsilon = 0.0;
    values = NULL;
    oldValues = NULL;
    values_tn_1 = NULL;
    values_tn_2 = NULL;
    values_tn_3 = NULL;
    regions.Clear();
    nodes.Clear();
    elements.Clear();
    neighInputRegions.Clear();
    materials.Clear();
  }

  PDECoupling::CouplingInterface &
  PDECoupling::CouplingInterface::operator= (const CouplingInterface & x)
  {
    regions = x.regions;
    regionType = x.regionType;
    nodes = x.nodes;
    elements = x.elements;
    neighInputRegions = x.neighInputRegions;
    materials = x.materials;

    if (x.values != NULL) {
      if (x.values->GetEntryType() == BaseMatrix::COMPLEX) {
        std::cerr << "IM FALSCHEN TEIL" << std::endl;
        EXCEPTION( "Sorry, wrong branch :) ");
        values = new Vector<Complex>(dynamic_cast<Vector<Complex>&>(*(x.values)));
        oldValues = new Vector<Complex>(dynamic_cast<Vector<Complex>&>(*(x.oldValues)));
        values_tn_1 = new Vector<Complex>(dynamic_cast<Vector<Complex>&>(*(x.values_tn_1)));
        values_tn_2 = new Vector<Complex>(dynamic_cast<Vector<Complex>&>(*(x.values_tn_2)));
        values_tn_3 = new Vector<Complex>(dynamic_cast<Vector<Complex>&>(*(x.values_tn_3)));

      } else {
        values = new Vector<Double>(dynamic_cast<Vector<Double>&>(*(x.values)));
        oldValues = new Vector<Double>(dynamic_cast<Vector<Double>&>(*(x.oldValues)));
        values_tn_1 = new Vector<Double>(dynamic_cast<Vector<Double>&>(*(x.values_tn_1)));
        values_tn_2 = new Vector<Double>(dynamic_cast<Vector<Double>&>(*(x.values_tn_2)));
        values_tn_3 = new Vector<Double>(dynamic_cast<Vector<Double>&>(*(x.values_tn_3)));
      }
    }  

    dof = x.dof;
    numNodes = x.numNodes;
    numElems = x.numElems;
    normtype = x.normtype;
    epsilon = x.epsilon;
    
    return *this;

  }
  
  
  template<class Archive>
  void PDECoupling::CouplingInterface::serialize(Archive & ar, 
                                                 const unsigned int version) {
    
    EXCEPTION( "Not implemented at the moment!" );

    // The problem at this point is that between reading and writing a 
    // restart file we can not keep the pointers contained in 
    // elements. 
    // Anyway, the whole mechanism has to be rewritten when we incorporate
    // the concept of p-fem, where we do not have only values associated
    // with nodes, but also edges surfaces and so on.
    // In this case we do have even more problems, as we cannot repsresent
    // these quantities without having the eqnMap object of the related
    // PDE available.
  }

  PDECoupling::CouplingInterface::
  CouplingInterface(const CouplingInterface &x ){
    *this = x;
 
  }

  PDECoupling::CouplingInterface::~CouplingInterface()
  {
    if(values != NULL)
      delete values;
    if(oldValues != NULL)
      delete oldValues;
    if(values_tn_1 != NULL)
      delete values_tn_1;
    if(values_tn_2 != NULL)
      delete values_tn_2;
    if(values_tn_3 != NULL)
      delete values_tn_3;
  }


  PDECoupling::PDECoupling(Grid * aptgrid)
  {
  
      ptGrid_ = aptgrid;
  
    defaultEpsilon = 1e-5;
    defaultNormType = L2REL;
  
  }

  PDECoupling::~PDECoupling()
  {
  
    // NOTE: since each interface exists only one time
    // but is referenced two times (once as input and once
    // as output), it must be deleted only once. In this 
    // case the input-interface pointer was used.
  
    //for (UInt i=0; i<outputInterfaces_.GetSize(); i++)
    //if (outputInterfaces_[i]) 
    //  delete outputInterfaces_[i];
  
    for (UInt i=0; i<inputInterfaces_.GetSize(); i++)
      if (inputInterfaces_[i]) 
        delete inputInterfaces_[i];
  
  }

  void PDECoupling::RegisterInput(CouplingInputType InType, SolutionType Quantity)

  {
    
    inputTypes_.Push_back(InType);
    inputQuantities_.Push_back(Quantity);  
    inputInterfaces_.Push_back(NULL);
  }


  void PDECoupling::AddInput(SolutionType quantity, 
                             StdVector<std::string> &region, 
                             CouplingRegionType regionType,
                             StdVector<std::string> &neighRegions,
                             Double epsilon,
                             NormType normtype,
                             StdVector<PDECoupling*> & couplings)
  {
  
  
    // Logging output
    LOG_DBG(pdecpl) << "Adding new input";
    LOG_DBG(pdecpl) << "\tQuantity: " << SolutionTypeEnum.ToString(quantity);
    LOG_DBG(pdecpl) << "\tRegions: " << region.Serialize();
    LOG_DBG(pdecpl) << "\tNeighbor Regions: " << neighRegions.Serialize();
    LOG_DBG(pdecpl) << "";
    LOG_DBG(pdecpl) << "";
    LOG_DBG(pdecpl) << "";
    
    CouplingInterface *myInterface = 0; 
    Integer myNum = -1;
  
    // search matching  Quantity
    for (UInt i=0; i<inputQuantities_.GetSize(); i++)
      {
        if (inputQuantities_[i] == quantity)
          myNum = i; 
      }
  
    if (myNum == -1)
      {
        EXCEPTION( "Quantity '" << SolutionTypeEnum.ToString(quantity) 
                   <<  "' not registered for PDE '"
                   << myPDE_->GetName() << "'" );
      }
  
  
    // Add coupling input as output to pde, which can calculate specified quantity
    CouplingOutputType myOutputType = Input2OutputType(inputTypes_[myNum]);
  
    UInt i=0;
    while (myInterface == 0 && i<couplings.GetSize())
      {
        myInterface = couplings[i++]->AddOutput(myOutputType, quantity, 
                                                region, regionType);
      }
  
    // If no pde has the specified quantity as output
    if (myInterface == NULL)
      {
        EXCEPTION( "Qantity \'" << SolutionTypeEnum.ToString(quantity) 
                   << "\' can not be calculated with current set of PDEs" );
      }
  
    // i now contains the number of
    // the coupling object, which shares the data
    // of this coupling object
    i--;  // because i has been incremented before;
    StdPDE * oppositePDE = couplings[i]->myPDE_;
  
  
    // Set dof according to myPDE's dof
    // if the inputType is COORD (= moving geometry)
    // then the number of dofs is equal to the
    // dimension of the grid, otherwise the number dofs
    // must be equal to that of my PDE
    if (inputTypes_[myNum] == COORD || 
        inputTypes_[myNum] == ID_BC || 
        inputTypes_[myNum] == GRID_VEL)
      myInterface->dof = myPDE_->dim_;
    else
      myInterface->dof = myPDE_->results_[0]->dofNames.GetSize();
  
    // Initialize the values and oldValues arrays
    //myInterface->values->SetDof(myPDE_->dofspernode_);
    //myInterface->values->Init();
    //myInterface->oldValues->SetDof(myPDE_->dofspernode_);
    //myInterface->oldValues->Init();
  
    // set normtype and epsilon
    //myInterface->epsilon = defaultEpsilon;
    //conf->ifget(inputQuantities_[myNum], myInterface->epsilon, "coupling", "tolerance"); 

    myInterface->epsilon = epsilon;
    myInterface->normtype = normtype;

    // set neighbouring region name(s)
    // if only entry is "all", then add all regions
    UInt numAllOccured = 0;

    if (neighRegions.GetSize() != 0){
      for (UInt iRegion=0; iRegion<neighRegions.GetSize(); iRegion++) {
        if (neighRegions[iRegion] == "all")
          numAllOccured ++;
      }
      if (numAllOccured == neighRegions.GetSize())
        {
          neighRegions.Clear();
          // get paramnode of oppsite pde
          ParamNodeList regionNodes = 
            oppositePDE->GetParamNode()->Get("regionList")->GetList("region");
          for( UInt iRegion = 0; iRegion < regionNodes.GetSize(); iRegion++ ) {
            neighRegions.Push_back( regionNodes[iRegion]->Get("name")->As<std::string>() );
          }
        }
      else if (numAllOccured > 0)
        {
          EXCEPTION( "AddInputCoupling: For 'neighbourRegion' either only 'all' "
                     << "or only region names are allowed!" );
        }
    
    }
    myInterface->neighInputRegions = neighRegions;
  
    if (myInterface->elements.GetSize() != 0
        && myInterface->regionType == SURFACE)
      {
        // Set the material of the interface.
        // Since the Inputcoupling values are computed by another PDE
        // as Output values, it might need the material paramter
        // for some integrators.
      
        // 1. Step: get the neighbouring elements 
      
        StdVector<Elem*> & interfaceElems = myInterface->elements;
        StdVector<Elem*>  actSubdomain;
        StdVector<Elem*> possibleNeighbours, neighbours;
    
      
      
        // for (UInt iSd=0; iSd < myPDE_->subdoms_.GetSize(); iSd++)
        //      {GetVolNeighboursForSurf
        //        ptGrid_->GetElemSD(actSubdomain, myPDE_->subdoms_[iSd]);
        //        for (UInt j=0; j<actSubdomain.GetSize(); j++)
        //          possibleNeighbours.Push_back(actSubdomain[j]);
        //      }
      
        //  ptGrid_->GetElemsNextToSurface( myInterface->neighbours, *interfaceElems,
        //                                    myPDE_->subdoms_ );

        //       if (!myInterface->neighbours.GetSize())
        //      EXCEPTION("No neighbours for element coupling found!");
      
      
      
      
        //       // 2. Step: For each interface element, set the
        //       //          material parameter according to its neighbour
      
        //       for (UInt i=0; i<myInterface->neighbours.GetSize(); i++)
        //      {
        //        bool subdomFound = false;
        //        UInt subDomNr = 0;
          
        //        for (subDomNr=0; subDomNr<myPDE_->subdoms_.GetSize(); subDomNr++)
        //          if (myPDE_->subdoms_[subDomNr] == myInterface->neighbours[i]->regionId)
        //            {
        //              subdomFound = true;
        //              break;
        //            }
          
          
        //        if (subdomFound == false)
        //          EXCEPTION("Subdomain name of neighbouring elements was not found");
          
        //        myInterface->materials[i] = &(myPDE_->materialData_[subDomNr]);           
        //      }
      
      
      
        // Set the material of the neighbours elements of the "opposite" PDE
        // 1. Step: get the neighbouring elements
        //       possibleNeighbours.Clear();
      
      
        //       for (UInt iSd=0; iSd < oppositePDE->subdoms_.GetSize(); iSd++)
        //         {
        //           ptGrid_->GetVolElems(actSubdomain, oppositePDE->subdoms_[iSd]);
        //        for (UInt j=0; j<actSubdomain.GetSize(); j++)
        //          possibleNeighbours.Push_back(actSubdomain[j]);
        //      }
      
        //       ptGrid_->GetElemsNextToSurface( neighbours,
        //                                       *interfaceElems,
        //                                       oppositePDE->subdoms_ );
      
        //        if (!neighbours.GetSize())
        //      EXCEPTION("No opposite neighbours for element coupling found!");
      
        //       for (UInt i=0; i<neighbours.GetSize(); i++)
        //         {
        //           bool subdomFound = false;
        //           UInt subDomNr = 0;
          
        //           for (subDomNr=0; subDomNr<oppositePDE->subdoms_.GetSize(); subDomNr++)
        //             if (oppositePDE->subdoms_[subDomNr] == neighbours[i]->regionId)
        //               {
        //                 subdomFound = true;
        //                 break;
        //               }
                  
          
        //        if (subdomFound == false)
        //          EXCEPTION("Subdomain name of neighbouring elements was not found");
          
        //        myInterface->oppositePdeMaterials[i] = &(oppositePDE->materialData_[subDomNr]);
        //      }
      
        Integer index = -1;
        myInterface->oppositePdeMaterials.Resize(interfaceElems.GetSize());
        myInterface->oppositePdeMaterials.Init(NULL);

        for ( UInt iElem = 0; iElem < interfaceElems.GetSize(); iElem++ ) {
        
          SurfElem const & myElem = 
            dynamic_cast<SurfElem &>(*interfaceElems[iElem]);

          index = myPDE_->subdoms_.Find(myElem.ptVolElem1->regionId);

          if ( index == -1 && (myElem.ptVolElem2!=NULL)) {
            index = myPDE_->subdoms_.Find(myElem.ptVolElem2->regionId);
          }
        
          if ( index == -1 && (myElem.ptVolElem2!=NULL)) {
            EXCEPTION( "PDECoupling::AddOutput: For Surface element Nr. " 
                       << " I found no according region in PDE '"
                       << myPDE_->GetName() << "'!" );
          }
          //In case we have only one volume elem neighbor we assume same mat index
          if ( index == -1 && (myElem.ptVolElem2==NULL)) {
          myInterface->oppositePdeMaterials[iElem] = 
            (myPDE_->materials_.begin()->second);
          }
          else
            myInterface->oppositePdeMaterials[iElem] = 
              (myPDE_->materials_[myPDE_->subdoms_[index]]);
        
        }
      
      } // end if
  
    inputInterfaces_[myNum] = myInterface;
  
  
  }
  
  void PDECoupling::Finalize() {
    
    // Create copy of output interfaces
    StdVector<CouplingOutputType> tempOutputTypes;     
    StdVector<SolutionType> tempOutputQuantities;
    StdVector<CouplingInterface*> tempOutputInterfaces;
    for( UInt i = 0; i < outputInterfaces_.GetSize(); ++i ) {
      if( outputInterfaces_[i] != NULL) {
        tempOutputTypes.Push_back(outputTypes_[i]);
        tempOutputQuantities.Push_back(outputQuantities_[i]);
        tempOutputInterfaces.Push_back(outputInterfaces_[i]);
      }
    }
    // store arrays back
    outputTypes_ = tempOutputTypes;
    outputQuantities_ = tempOutputQuantities;
    outputInterfaces_ = tempOutputInterfaces;
      
    // Create copy of input interfaces
    StdVector<CouplingInputType> tempInputTypes;     
    StdVector<SolutionType> tempInputQuantities;
    StdVector<CouplingInterface*> tempInputInterfaces;
    for( UInt i = 0; i < inputInterfaces_.GetSize(); ++i ) {
      if( inputInterfaces_[i] != NULL) {
        tempInputTypes.Push_back(inputTypes_[i]);
        tempInputQuantities.Push_back(inputQuantities_[i]);
        tempInputInterfaces.Push_back(inputInterfaces_[i]);
      }
    }
    // store array back
    inputTypes_ = tempInputTypes;
    inputQuantities_ = tempInputQuantities;
    inputInterfaces_ = tempInputInterfaces;
    
  }

  PDECoupling::CouplingInterface* PDECoupling::AddOutput(CouplingOutputType outputType, 
                                                         SolutionType quantity, 
                                                         StdVector<std::string> & regions,
                                                         CouplingRegionType regionType)
  {
  
    // search if output exists already
    for (UInt i=0; i<outputTypes_.GetSize(); i++)
      if (outputTypes_[i] == outputType 
          && outputQuantities_[i] == quantity
          && outputInterfaces_[i]->regions == regions)
        return outputInterfaces_[i];
  
  
    // if not find out if PDE can calculate desired output
    if (myPDE_->HasOutput(quantity) == false)
      return 0;
  
    // create new Coupling Output
    outputTypes_.Push_back(outputType);
    outputQuantities_.Push_back(quantity);

    CouplingInterface *myInterface = new CouplingInterface;  
    myInterface->regions = regions;
    myInterface->regionType = regionType;

    StdVector<UInt> nodesConverted;
    StdVector<Elem*> SD, tempSD;
    StdVector<SurfElem*> surfElems;
    StdVector<UInt> nodes;
    UInt inode = 0;
    UInt numNodes = 0;
    StdVector<Elem*> auxElems;
    std::string errMsg;
    StdVector<RegionIdType> regionIdVec;
  
    // Get Elements/nodes of coupling region
    SD.Clear();
    switch (regionType)
      {
      case REGION:

        ptGrid_->GetRegion().Parse(regions, regionIdVec);

        if( outputType == NODE ) {
          numNodes = ptGrid_->GetNumNodes( regionIdVec );
          myInterface->nodes.Reserve(numNodes);
          
          
          for (UInt iSD=0; iSD<regions.GetSize(); iSD++)
            {
              ptGrid_->GetNodesByRegion(nodes, regionIdVec[iSD]);
              for (UInt iNode=0; iNode<nodes.GetSize(); iNode++) {
                if( myInterface->nodes.Find(nodes[iNode]) == -1 )
                  myInterface->nodes.Push_back(nodes[iNode]);
              }
            }
          
          myInterface->numNodes = myInterface->nodes.GetSize();
          
          // Check if any nodes at all were found
          if ( myInterface->numNodes == 0){
            errMsg  = "Coupling::AddOutput: The region(s) ";
            for (UInt k=0; k<regions.GetSize()-1; k++) {
              errMsg += "'";
              errMsg += regions[k];
              errMsg += "', ";
            }
            errMsg +="'";
            errMsg += regions[regions.GetSize()-1];
            errMsg += "' were not found in the mesh or contain no nodes.\n";
            errMsg += "Please check you mesh file!";
            EXCEPTION(errMsg.c_str() );
          }
      
          //myInterface->values.resize(myInterface->nodes.GetSize());
          //myInterface->oldValues.resize(myInterface->nodes.GetSize());
          //  myInterface->values->SetNumNodes(myInterface->nodes.GetSize());
          //       myInterface->oldValues->SetNumNodes(myInterface->nodes.GetSize());
        } else {
          
          tempSD.Clear();
          for (UInt iSD=0; iSD<regions.GetSize(); iSD++) {
            ptGrid_->GetElems( tempSD, regionIdVec[iSD] );

            for (UInt iElem=0; iElem<tempSD.GetSize(); iElem++) {
              myInterface->elements.Push_back( tempSD[iElem] );
              }
            }
          
          myInterface->numElems = myInterface->elements.GetSize();
        
        }
        
        break;

      case NODES:
      

        // count complete number of nodes
        for (UInt i=0; i<regions.GetSize(); i++)
          numNodes+= ptGrid_->GetNumNodes(regions[i]);

        myInterface->nodes.Resize(numNodes);
        myInterface->nodes.Init();
      
        inode =0;

        // get for each nodeslist all nodes
        for (UInt i=0; i<regions.GetSize(); i++)
          {
            ptGrid_->GetNodesByName(nodesConverted, regions[i]);
          
            for ( UInt k=0; k< nodesConverted.GetSize(); k++, inode++)
              myInterface->nodes[inode] = nodesConverted[k];
          }

        myInterface->numNodes = myInterface->nodes.GetSize();

        // Check if any nodes at all were found
        if ( myInterface->numNodes == 0){
          errMsg  = "Coupling::AddOutput: The nodes(s) ";
          for (UInt k=0; k<regions.GetSize()-1; k++) {
            errMsg += "'";
            errMsg += regions[k];
            errMsg += "', ";
          }
          errMsg +="'";
          errMsg += regions[regions.GetSize()-1];
          errMsg += "' were not found in the mesh or contain no nodes.\n";
          errMsg += "Please check you mesh file!";
          EXCEPTION(errMsg.c_str() );
        }
        break;

      case SURFACE:

        ptGrid_->GetRegion().Parse( regions, regionIdVec );
        numNodes = ptGrid_->GetNumNodes( regionIdVec );
        myInterface->nodes.Reserve(numNodes);

        for (UInt iSD=0; iSD<regions.GetSize(); iSD++) {

          // first get the surface elements
          ptGrid_->GetSurfElems( surfElems, regionIdVec[iSD] );
          for (UInt iElem=0; iElem<surfElems.GetSize(); iElem++)
            SD.Push_back(surfElems[iElem]);
        
          // then get the according coupling nodes
          ptGrid_->GetNodesByRegion(nodes, regionIdVec[iSD]);
        
          for (UInt iNode=0; iNode<nodes.GetSize(); iNode++) {
            myInterface->nodes.Push_back(nodes[iNode]);
          }
        
        }
      
      
        // Get the nodes from BCs
        //myInterface->elements = ptBCs_->getEdgesBC(region);
        myInterface->elements = SD;
      

        myInterface->numNodes = myInterface->nodes.GetSize();

        // Check if any nodes at all were found
        if ( myInterface->numNodes == 0){
          errMsg  = "Coupling::AddOutput: The surface(s) ";
          for (UInt k=0; k<regions.GetSize()-1; k++) {
            errMsg += "'";
            errMsg += regions[k];
            errMsg += "', ";
          }
          errMsg +="'";
          errMsg += regions[regions.GetSize()-1];
          errMsg += "' were not found in the mesh or contain no nodes.\n";
          errMsg += "Please check you mesh file!";
          EXCEPTION(errMsg.c_str() );
        }
      
        myInterface->numElems = myInterface->elements.GetSize();
        myInterface->materials.Resize(myInterface->elements.GetSize());
        myInterface->materials.Init( NULL );
            

        break;
      }

    outputInterfaces_.Push_back(myInterface);
  
    return myInterface;
  }




  void PDECoupling::SetPDE(StdPDE * aPDE)
  {
    myPDE_ = aPDE;
  }

  void PDECoupling::SetOutputNumNodes(UInt i, UInt nnodes)
  {
  
    outputInterfaces_[i]->numNodes = nnodes;
    
  }  

  void PDECoupling::SetOutputNumElems(UInt i, UInt nelems)
  {
    outputInterfaces_[i]->numElems = nelems;
  
  }  


  void PDECoupling::SetOutputDof(UInt i, UInt dof)
  {
    outputInterfaces_[i]->dof = dof;
  }

  StdPDE*  PDECoupling::GetPDE()
  {
    return myPDE_;
  }


  UInt PDECoupling::GetNumInputCouplings()
  {
    
    UInt numCouplings = 0;
    
    for ( UInt i = 0; i < inputInterfaces_.GetSize(); i++ ) {
      if ( inputInterfaces_[i] != NULL )
        numCouplings++;
    }
    
    return numCouplings;
  
  }


  UInt PDECoupling::GetNumOutputCouplings()
  {
    UInt numCouplings = 0;
    
    for ( UInt i = 0; i < outputInterfaces_.GetSize(); i++ ) {
      if ( outputInterfaces_[i] != NULL )
        numCouplings++;
    }
    return numCouplings;
//     return outputQuantities_.GetSize();
  }


  void PDECoupling::CreateCouplingVector(UInt i,
                                         bool isComplex)
  {

    UInt numNodes = outputInterfaces_[i]->numNodes;
    UInt numElems = outputInterfaces_[i]->numElems;
    UInt dof =  outputInterfaces_[i]->dof;

    if (outputInterfaces_[i]->values != NULL || 
        outputInterfaces_[i]->oldValues != NULL ||
        outputInterfaces_[i]->values_tn_1 != NULL ||
        outputInterfaces_[i]->values_tn_2 != NULL ||
        outputInterfaces_[i]->values_tn_3 != NULL)
      EXCEPTION("PDECoupling::CreateCouplingVector: This function may only\
           be called for Output-Coupling interfaces!" );

    if (isComplex)
      {
        outputInterfaces_[i]->values = new Vector<Complex>;
        outputInterfaces_[i]->oldValues = new Vector<Complex>;
      outputInterfaces_[i]->values_tn_1 = new Vector<Complex>;
      outputInterfaces_[i]->values_tn_2 = new Vector<Complex>;
      outputInterfaces_[i]->values_tn_3 = new Vector<Complex>;
      } 
    else
      {
        outputInterfaces_[i]->values = new Vector<Double>;
        outputInterfaces_[i]->oldValues = new Vector<Double>;
      outputInterfaces_[i]->values_tn_1 = new Vector<Double>;
      outputInterfaces_[i]->values_tn_2 = new Vector<Double>;
      outputInterfaces_[i]->values_tn_3 = new Vector<Double>;
      }
  
    // Determine size of vector
    UInt size = 0;
    
    if ( numElems != 0 ) {
      size = numElems;
    }
    if ( numNodes != 0 ) {
      size = numNodes;
    }
    if ( size == 0 ) {
      EXCEPTION( "Will not allocate a vector of size 0" );
    }

    // resize vector with coupling values
    outputInterfaces_[i]->values->Resize(dof*size);
    outputInterfaces_[i]->values->Init( );
    outputInterfaces_[i]->oldValues->Resize(dof*size);
    outputInterfaces_[i]->oldValues->Init( );
    outputInterfaces_[i]->values_tn_1->Resize(dof*size);
    outputInterfaces_[i]->values_tn_1->Init();
    outputInterfaces_[i]->values_tn_2->Resize(dof*size);
    outputInterfaces_[i]->values_tn_2->Init();
    outputInterfaces_[i]->values_tn_3->Resize(dof*size);
    outputInterfaces_[i]->values_tn_3->Init();

  }


  CouplingOutputType PDECoupling::Input2OutputType(CouplingInputType inputType)
  {

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

      case GRID_VEL:
        // Grid Velocity coupling -> node values needed
        return NODE;
        break;
      }

    EXCEPTION( "No output found for input" );
    
    // 

  }

  void PDECoupling::GetMemento(CouplingMemento & memento) {
  
    memento.inputTypes_ = inputTypes_;
    memento.inputQuantities_ = inputQuantities_;
    memento.inputInterfaces_.Resize(inputInterfaces_.GetSize());

    for (UInt i=0; i<inputInterfaces_.GetSize(); i++){
      memento.inputInterfaces_[i] = *(inputInterfaces_[i]);
    }
    memento.isSet_ = true;
  }


  void PDECoupling::SetMemento(CouplingMemento & memento) {
  
    std::string errMsg, warnMsg;

    // if there is no information in the memento just leave
    if (memento.isSet_ == false)
      return;
  
    // Check if size of internal vectors
    if ((memento.inputTypes_.GetSize() != memento.inputQuantities_.GetSize()) ||
        (memento.inputTypes_.GetSize() != memento.inputInterfaces_.GetSize())) {
      errMsg  = "PDECoupling::SetMemento: Incnsistent definition!\n";
      errMsg += "Size of 'inputTypes_': ";
      errMsg += inputTypes_.GetSize();
      errMsg += "\nSize of 'inputQuantities_': ";
      errMsg += inputQuantities_.GetSize();
      errMsg += "\nSize of 'inputInterfaces_: ";
      errMsg += inputInterfaces_.GetSize();
      EXCEPTION(errMsg.c_str() );
    }

    // Iterate over all quantities and check which values can be assigned

    // loop over all memento interfaces
    for (UInt iMem=0; iMem<memento.inputTypes_.GetSize(); iMem++) {
      bool couplingFound = false;

      // loop over all own interfaces
      for (UInt iOwn=0; iOwn<inputTypes_.GetSize(); iOwn++) {
        if ((memento.inputTypes_[iMem] == inputTypes_[iOwn]) && 
          
            (memento.inputQuantities_[iMem] == inputQuantities_[iOwn]) &&
          
            (memento.inputInterfaces_[iMem].regions  == 
             inputInterfaces_[iOwn]->regions) &&
          
            (memento.inputInterfaces_[iMem].regionType  == 
             inputInterfaces_[iOwn]->regionType) &&
          
            (memento.inputInterfaces_[iMem].nodes  == 
             inputInterfaces_[iOwn]->nodes) &&
          
            (memento.inputInterfaces_[iMem].elements  == 
             inputInterfaces_[iOwn]->elements) &&
          
            (memento.inputInterfaces_[iMem].dof  == 
             inputInterfaces_[iOwn]->dof) &&
          
            (memento.inputInterfaces_[iMem].numNodes  == 
             inputInterfaces_[iOwn]->numNodes) &&
          
            (memento.inputInterfaces_[iMem].numElems  == 
             inputInterfaces_[iOwn]->numElems) &&
          
            (memento.inputInterfaces_[iMem].values->GetEntryType() == 
             inputInterfaces_[iOwn]->values->GetEntryType()))
          {
          
            // This section is reached, if the two couplings
            // are the same
            if (inputInterfaces_[iOwn]->values->GetEntryType()
                == BaseMatrix::COMPLEX){
            
              // --- Complex values --
              dynamic_cast<Vector<Complex>&> (*(inputInterfaces_[iOwn]->values)) 
                = dynamic_cast<Vector<Complex>&>
                (*(memento.inputInterfaces_[iMem].values));
            
            } else {
              // --- Real values --

              dynamic_cast<Vector<Double>&> (*(inputInterfaces_[iOwn]->values))
                = dynamic_cast<Vector<Double>&>
                (*(memento.inputInterfaces_[iMem].values));
            }
            
            couplingFound = true;
            break;
          } // if statement
      
        if (! couplingFound) {
          warnMsg  = "PDECoupling::SetMemento: The coupling quantitiy '";
          warnMsg += SolutionTypeEnum.ToString(memento.inputQuantities_[iMem]);
          warnMsg += "' at the interface(s) '";
          for (UInt i=0; i<memento.inputInterfaces_[iMem].regions.GetSize()-1; i++) {
            warnMsg +=  memento.inputInterfaces_[iMem].regions[i];
            warnMsg += ", ";
          }
          warnMsg += memento.inputInterfaces_[iMem].
            regions[memento.inputInterfaces_[iMem].regions.GetSize()-1];
          warnMsg +="' is not present anymore in the current coupling set for PDE '";
          warnMsg += myPDE_->GetName();
          warnMsg += "'.";
          WARN(warnMsg);
        } //if coupling not found
      } // loop over own couplings
    } // looop over memento couplings
  
  }

  //   std::ostream& operator << ( std::ostream & out , const CouplingInterface & inter)
  //   {
  //     out << myendl 
  //      << "Interface: coupling region: " << inter.region << myendl
  //      << "           nr. coupling nodes: " <<  inter.nodes.GetSize() << myendl
  //      << "           nr. coupling elements: " <<  inter.elements.GetSize() << myendl
  //      << "           name of 1. coupling element: " <<  inter.elements[0].namesd << myendl;
  //      << "           nr. neighbours elements: " <<  inter.neighbours.GetSize() << myendl
  //      << "           name of 1. neighbours element: " <<  inter.neighbours[0].namesd << myendl;

  //   } 

} // end of namespace

#include "boost/serialization/export.hpp"

BOOST_CLASS_EXPORT_GUID(CoupledField::PDECoupling::CouplingInterface,
                        "CoupledField_CouplingInterface")
