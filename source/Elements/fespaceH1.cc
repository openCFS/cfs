#include "fespaceH1.hh"
#include "H1Elems.hh"

namespace CoupledField {

  //! Constructor
  FeSpaceH1::FeSpaceH1() : FeSpace() {
    type_ = H1_LO;
  }
  
  FeSpaceH1::~FeSpaceH1(){
  }

  //! Set type of basis
  void FeSpaceH1::SetBasis( BasisType type ){
    basis_ = type;
  }
  
  void FeSpaceH1::AddFeFunction( shared_ptr<BaseFeFunction> fct ){
    feFunction_ = fct;
    return;
  }
  
  BaseFE* FeSpaceH1::GetFe( const EntityIterator ent ){
    if(refElems_.find(ent.GetElem()->type) == refElems_.end()){
      EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
    }
    return refElems_[ent.GetElem()->type];
  }

  //! Returns the number of (vectorial) unkowns on the element
  UInt FeSpaceH1::GetNumFunctions( const EntityIterator ent ){
    //just for debugging purpose
    if(refElems_.find(ent.GetElem()->type) == refElems_.end()){
      EXCEPTION("FeSpaceH1::GetNumFunctions(const EnitityIterator): requested fetype which is noch supported by space");
    }
    //return refElems_[ent.GetElem()->ptElem->feType()]->GetNumFncs(feFunction_->GetResultInfo()->fctType);
    return refElems_[ent.GetElem()->type]->GetNumFncs();
  }

  //! Return equation numbers for a all DOFs
  void FeSpaceH1::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ){

    //Get result for the feFunction
    shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();

    //Get "dimension" of one Unknown
    UInt dofsPerUnknown = feFctResult->dofNames.GetSize();

    //This if clause should be avoided if the functionality for higher order entities
    //is implemented
    if( feFctResult->definedOn == ResultInfo::NODE ) {

      //First cover the nodal Case
      if ( ent.GetType() == EntityList::NODE_LIST ) {
        UInt node = ent.GetNode();
        eqns.Resize(dofsPerUnknown);
        eqns.Init();
        for(UInt iDof = 0; iDof < dofsPerUnknown; iDof++){
          eqns[iDof] = nodeMap_[node][iDof];
        }
      } else if( ent.GetType() == EntityList::ELEM_LIST ||
                 ent.GetType() == EntityList::SURF_ELEM_LIST){
        StdVector<UInt> nodes;
        GetNodesOfElement(nodes,ent.GetElem());
        eqns.Resize( nodes.GetSize() * dofsPerUnknown );
        eqns.Init();
        for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
          for(UInt iDof = 0; iDof < dofsPerUnknown; iDof++){
            eqns[(iNode*dofsPerUnknown) + iDof] = nodeMap_[nodes[iNode]][iDof];
          }
        }
      } else {
        EXCEPTION("In FeSpaceH1::GetEqns():  Supplied an iterator which is not supported by FeSpace");
      }
    }else{
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //TODO: Higher order Dofs (Edges, Faces, Interior)
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      EXCEPTION("In FeSpaceH1::GetEqns(): This FeSpace currently only supports nodal equation numbering");
    }
  }

  //! Return equation numbers for a specific dof
  void FeSpaceH1::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                        , UInt dof ){ 
    //Get result for the feFunction
    shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();

    //Get "dimension" of one Unknown
    UInt dofsPerUnknown = feFctResult->dofNames.GetSize();

    //First cover the nodal/grid case
    if ( ent.GetType() == EntityList::NODE_LIST ) {
      UInt node = ent.GetNode();
      eqns.Resize(1);
      eqns.Init();
      eqns[0] = nodeMap_[node][dof];

    } else if( ent.GetType() == EntityList::ELEM_LIST ||
               ent.GetType() == EntityList::SURF_ELEM_LIST){
      StdVector<UInt> nodes;
      GetNodesOfElement(nodes,ent.GetElem());
      eqns.Resize( nodes.GetSize() );
      eqns.Init();
      for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
        eqns[iNode] = nodeMap_[nodes[iNode]][dof];
      }
    } else {
      EXCEPTION("In FeSpaceH1::GetEqns(StdVector,EntityIterator,UInt):  Supplied an iterator which is not supported by FeSpace");
    }
  
  }

  //! get a nodal equation number
  UInt FeSpaceH1::GetNodeEqn(UInt nodeNr, UInt dof){
    return nodeMap_[nodeNr][dof];
  }

  //! Get Equation numbers for a specific element
  void FeSpaceH1::GetElemEqns(StdVector<Integer>& eqns,const Elem* elem){
  }

  //! Get Equation numbers for a specific element
  void FeSpaceH1::GetElemEqns(StdVector<Integer>& eqns,const Elem* elem, UInt dof){
  }

  //! Map equations i.e. intialize object
  void FeSpaceH1::Finalize(){

    UInt numEqns_ = 0;
    UInt numFreeEquations_ = 0;
    CreateVirtualNodes();
    //Determine boundary Unknowns
    MapBCs();
    MapNodalEqns(1);
    MapNodalEqns(2);
    isFinalized_ = true;
  }

  //! Map BC Equation 
  void FeSpaceH1::MapBCs(){

  
  }
  
  //=======================================================
  //Perform Nodal equation Numbering
  //=======================================================
  void FeSpaceH1::MapNodalEqns(UInt phase){
    UInt dofsPerUnknown = 0;
    shared_ptr<ResultInfo> feFctResult;
    StdVector< shared_ptr<EntityList> > fctEntList;
    StdVector< shared_ptr<EntityList> >::iterator entIt;
    std::map< Integer , StdVector<BC_Type> >::iterator bcIt;
    UInt actNode = -1;

    switch(phase){
      case 1:
        // number the nodal equations if they are not contained in the BcKeys
        // due to our virtual node array we have a valid structure

        //Get result for the feFunction
        feFctResult = feFunction_->GetResultInfo();
        // Get number of dofs
        dofsPerUnknown = feFctResult->dofNames.GetSize();
        
        //now loop over all nodes and assign an equation number
        for(UInt curNode = 0; curNode < nodes_.GetSize(); curNode++){
          actNode = nodes_[curNode];

          if(nodeMap_.eqns.find(actNode) == nodeMap_.eqns.end()){
            nodeMap_[actNode] = StdVector<Integer>(dofsPerUnknown,-1);
          }
          numUnknowns_++;

          for(UInt iDof = 0; iDof < dofsPerUnknown; iDof ++){
            if(nodeMap_.BcKeys.find(actNode) != nodeMap_.BcKeys.end()){
              if(nodeMap_.BcKeys[actNode][iDof] == NOBC){
                nodeMap_[actNode][iDof] = ++numEqns_;
                numFreeEquations_++;
              }
            }else if(nodeMap_[actNode][iDof] == -1){
              nodeMap_[actNode][iDof] = ++numEqns_;
              numFreeEquations_++;
            }
          }
        }
        break;
      case 2:
        bcIt = nodeMap_.BcKeys.begin();
        while(bcIt != nodeMap_.BcKeys.end()){
          //TODO: take the kind of BC into account
          //OPTIONAL: Delete the BcKeys Map afterwards?!
          for(UInt iDof =0; iDof < bcIt->second.GetSize(); iDof++){
            if(bcIt->second[iDof]!=NOBC){
               //nodeMap_[bcIt->first] = ++numEqns_;
              if(bcIt->second[iDof] == IDBC){
                nodeMap_[bcIt->first] = ++numEqns_;
              } else if(bcIt->second[iDof] == HDBC){
                nodeMap_[bcIt->first] = 0 ;
              } 
              //else if(bcIt->second[iDof] == CONSTRAINT){
              //  Integer masterNode = nodeMap_.slaveMasterNodes[bcIt->first];
              //  nodeMap_[bcIt->first] = -1*masterNode;
              //}
            }
          }
          bcIt++;
        }
        {
          std::map<std::pair<Integer,Integer>, std::pair<Integer,Integer>  >::iterator  conNodes = nodeMap_.constraintNodes.begin();
          for(conNodes; conNodes !=nodeMap_.constraintNodes.end() ; conNodes++){
            std::pair<Integer,Integer> slaveP = conNodes->first;
            std::pair<Integer,Integer> masterP = conNodes->second;
            nodeMap_[slaveP.first][slaveP.second] = -1*nodeMap_[masterP.first][masterP.second];
          }
        }
        break;
      default:
        EXCEPTION("In FeSpaceH1::MapNodalEqns(): Supplied wrong argument for the numbering phase");
        break;
    }
  }

  //! Reorder the equation Map (just for comptibility)
  void FeSpaceH1::ReorderEqnMap( StdVector<UInt> newOrder ){
    if(!newOrder.GetSize()){
      return;
    }
    //NODAL PART
    std::map< Integer , StdVector<Integer> >::iterator mapIt = nodeMap_.eqns.begin();
    while(mapIt != nodeMap_.eqns.end()){
      for(UInt iDof =0; iDof < mapIt->second.GetSize(); iDof++){
        if(mapIt->second[iDof] > 0){
          mapIt->second[iDof] = (Integer)newOrder[mapIt->second[iDof]-1];
        }else if (mapIt->second[iDof] < 0){
          mapIt->second[iDof] = -(Integer)newOrder[-mapIt->second[iDof]-1];
        }
      }
      mapIt++;
    }
  }

  void FeSpaceH1::PrintEqnMap(){
     std::map< Integer , StdVector<Integer> >::iterator mapIt = nodeMap_.eqns.begin();
     std::cout << "EQUATION MAPPING" << std::endl << std::endl;
     std::cout << "nodeNr \t|\t equationNr" << std::endl;
     std::cout << "-----------------------------------------------" << std::endl;
     while(mapIt != nodeMap_.eqns.end()){
       for(UInt iDof =0; iDof < mapIt->second.GetSize(); iDof++){
         std::cout << mapIt->first << "\t|\t" << mapIt->second[iDof] << std::endl;
       }
       mapIt++;
     }
  }


  //! Map Edge Equation Numbers
  void FeSpaceH1::MapEdgeEquations(UInt phase){
    EXCEPTION("In FeSpaceH1::MapEdgeEquations(UInt phase): This FeSpace currently only supports nodal equation numbering");
  }

  //! Map Face Equation Numbers
  void FeSpaceH1::MapFaceEquations(UInt phase){
    EXCEPTION("In FeSpaceH1::MapFaceEquations(UInt phase): This FeSpace currently only supports nodal equation numbering");
  }

  //! Map Inerior Equation Numbers
  void FeSpaceH1::MapInteriorEquations(UInt phase){
    EXCEPTION("In FeSpaceH1::MapInteriorEquations(UInt phase): This FeSpace currently only supports nodal equation numbering");
  }

  

}
