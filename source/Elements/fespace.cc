// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "fespace.hh"


namespace CoupledField {

  //! Constructor
  FeSpace::FeSpace(){
    
    isFinalized_ = false;
    numEqns_ = 0;
    numUnknowns_ = 0;
    numFreeEquations_ = 0;
    mapType_ = GRID;

    bcCounter_[NOBC] = 0;
    bcCounter_[HDBC] = 0;
    bcCounter_[IDBC] = 0;
    bcCounter_[CONSTRAINT] = 0;

    //now call the function defining its integration scheme

  }

  FeSpace::~FeSpace(){
  }

  void FeSpace::GetNodesOfEntities( StdVector<UInt>& nodes,
                                    shared_ptr<EntityList> ent,
                                    ElementEntityType entType ) {

    if (mapType_ == GRID){ // get name of entitylist
      std::string name= ent->GetName();
      feFunction_->GetGrid()->GetNodesByName( nodes, name );
    } else if (mapType_ == POLYNOMIAL) {
      EntityList::ListType defType = ent->GetType();
      EntityIterator entIt = ent->GetIterator();
      nodes.Resize(0);
      UInt nCount = 0; 
      switch (defType){
        case EntityList::REGION_LIST:
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            StdVector<Elem*> eList;
            feFunction_->GetGrid()->GetElems(eList,entIt.GetRegion());
            StdVector<UInt> eNodes;
            for (UInt i = 0; i < eList.GetSize(); i ++ ){

              GetNodesOfElement(eNodes,eList[i],entType);
              for (UInt aNode = 0; aNode < eNodes.GetSize(); aNode += 1 ){
                nodes.Push_back(eNodes[aNode]);
              }
            }
          }
        break;
        case EntityList::ELEM_LIST:
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            StdVector<UInt> eNodes;
            GetNodesOfElement(eNodes,entIt.GetElem(),entType);
            for (UInt aNode = 0; aNode < eNodes.GetSize(); aNode += 1 ){
              nodes.Push_back(eNodes[aNode]);
            }
          }
        break;
        case EntityList::SURF_ELEM_LIST:
          Warning(" FeSpaceH1::GetNodesOfEntities(): Going to treat a SURF_ELEM_LIST as a ELEM_LIST...");
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            StdVector<UInt> eNodes;
            GetNodesOfElement(eNodes,entIt.GetElem(),entType);
            for (UInt aNode = 0; aNode < eNodes.GetSize(); aNode += 1 ){
              nodes.Push_back(eNodes[aNode]);
            }
          }
        break;
        case EntityList::NODE_LIST:
          nCount = 0;
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            if(gridToVirtualNodes_.find(entIt.GetNode()) != gridToVirtualNodes_.end()){
              nodes.Push_back( gridToVirtualNodes_[entIt.GetNode()]);
            }
          }
        break;
        default:
          EXCEPTION("FeSpace::GetNodesOfEntities: Called with a list which is not supported ");
          break;
      }
    }
  }

  void FeSpace::GetNodesOfElement( StdVector<UInt>& nodes,
                          const Elem* ptElem,
                          ElementEntityType entType){
    UInt elemNum = ptElem->elemNum;
    if(virtualNodes_.find(elemNum) ==virtualNodes_.end()){
      EXCEPTION("FeSpace::GetNodesOfElement: Could not find requested element Number");
    }
    if(entType == ALL){
      nodes.Resize(virtualNodes_[elemNum][VERTEX].GetSize()+
                   virtualNodes_[elemNum][EDGE].GetSize()+
                   virtualNodes_[elemNum][FACE].GetSize()+
                   virtualNodes_[elemNum][INTERIOR].GetSize());
      ElemVirtualNodes::iterator nodeIt = virtualNodes_[elemNum].begin();
      UInt c = 0;
      while(nodeIt !=virtualNodes_[elemNum].end()){

        StdVector<UInt> entNodes =  nodeIt->second;
        for (UInt i = 0; i < entNodes.GetSize(); i++ ){
          nodes[c++] =  entNodes[i];
        }
        nodeIt++;
      }
    }else{
      nodes.Resize(virtualNodes_[elemNum][entType].GetSize());
      StdVector<UInt> entNodes =  virtualNodes_[elemNum][entType];
      for (UInt i = 0; i < entNodes.GetSize(); i++ ){
        nodes[i] =  entNodes[i];
      }
    }
  }

}
