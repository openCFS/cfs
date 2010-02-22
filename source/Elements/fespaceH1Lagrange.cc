// =====================================================================================
// 
//       Filename:  fespaceH1Lagrange.cc
// 
//    Description:  This implementes the space of lagrange elements the most
//                  important method here is the finalize method which fills the
//                  virtual node map.
// 
//        Version:  0.1
//        Created:  01/21/2010 10:42:37 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================
#include "Elements/fespaceH1Lagrange.hh"

namespace CoupledField{

    //! Constructor
    FeSpaceH1Lagrange::FeSpaceH1Lagrange(){
      mapType_ = GRID;
    }

    //! Destructor
    FeSpaceH1Lagrange::~FeSpaceH1Lagrange(){
    }

    void FeSpaceH1Lagrange::SetMapType( MappingType mapT){
      mapType_ = mapT;
      
      //build up the pointerMap
      
      // Note: at the moment this is not implemented in a 
      // clean fashion. The equation mapping currently relies
      // on the global node numbers of the element.
      // Instead, it should use the GetNumFncs(AnsatzFct::NODE)
      // method, which would deliver the correct number of
      // unknowns for the calculation element instead of the 
      // geometric element.
      if( mapType_ == GRID ) {
        refElems_[Elem::ET_LINE2]  = new FeH1LagrangeLine1();
        refElems_[Elem::ET_QUAD4]  = new FeH1LagrangeQuad1();
        refElems_[Elem::ET_HEXA8]  = new FeH1LagrangeHex1();
        refElems_[Elem::ET_LINE3]  = new FeH1LagrangeLine2();
        refElems_[Elem::ET_QUAD8]  = new FeH1LagrangeQuad2();
        refElems_[Elem::ET_HEXA20] = new FeH1LagrangeHex2();
      } else if (mapType_ == POLYNOMIAL) {
        refElems_[Elem::ET_LINE2]  = new FeH1LagrangeLineVar();
        refElems_[Elem::ET_QUAD4]  = new FeH1LagrangeQuadVar();
        refElems_[Elem::ET_HEXA8]  = new FeH1LagrangeHexVar();
        refElems_[Elem::ET_LINE3]  = new FeH1LagrangeLineVar();
        refElems_[Elem::ET_QUAD8]  = new FeH1LagrangeQuadVar();
        refElems_[Elem::ET_HEXA20] = new FeH1LagrangeHexVar();

        ////for lagrangian elements only isoparateric element orders are supported
        //std::map<FEType,BaseFE*>::iterator i = refElems_.begin();
        //for( ; i != refElems_.end(); ++i ) {
        //  i->second->SetIsoOrder(order_);
        //}

      }
    }

    //THIS IS COMMENTED OUT FOR NOW, BUR PERHAPS WE WILL NEED IT AGAIN WHEN THE GENERAL IMPLEMENTATION
    //IN THE FESPACEH1 CLASS IS NO LONGER SUFFICIENT
    ////! Return equation numbers
    //void FeSpaceH1Lagrange::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ){ 
    //}

    ////! Return equation numbers for a specific dof
    //void FeSpaceH1Lagrange::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
    //                      , UInt dof ){ 
    //}

    //! Map equations i.e. intialize object
    void FeSpaceH1Lagrange::Finalize(){
      /* Basic idea:
       * 1. Create the VirtualNode Array
       * 2. Map boundary conditions
       * 3. Map equations only based on the virtualNodeArray
      */

      UInt numEqns_ = 0;
      UInt numFreeEquations_ = 0;
      CreateVirtualNodes();
      //Determine boundary Unknowns
      MapBCs();
      MapNodalEqns(1);
      MapNodalEqns(2);
      isFinalized_ = true;
    }

    //! Reorder the equation Map (just for comptibility)
    void FeSpaceH1Lagrange::ReorderEqnMap( StdVector<UInt> newOrder ){
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

    //! print the equation map
    void FeSpaceH1Lagrange::PrintEqnMap(){
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

    //! Map BC Equation NUmbers
    void FeSpaceH1Lagrange::MapBCs(){
      StdVector<UInt> actNodes;
      //Get Grip of HdBC List for the fefunction
      const HdBcList hdbcs = feFunction_->GetHomDirichletBCs();
      HdBcList::const_iterator actHBC;

      UInt dofsPerUnknown = feFunction_->GetResultInfo()->dofNames.GetSize();

      for(actHBC = hdbcs.Begin(); actHBC != hdbcs.End(); actHBC++) {
        // Get EntityIterator
        GetNodesOfEntities(actNodes,(*actHBC)->entities);
        for(UInt iNode = 0 ; iNode < actNodes.GetSize();iNode++){
           if( nodeMap_.BcKeys.find(actNodes[iNode]) == nodeMap_.BcKeys.end()){
             //nodeMap_.BcKeys[node] = StdVector<BaseFeFunction::BC_Type>(dofsPerUnknown,BaseFeFunction::NOBC);
             nodeMap_.BcKeys[actNodes[iNode]] = StdVector<BC_Type>(dofsPerUnknown,NOBC);
           }
           nodeMap_.BcKeys[actNodes[iNode]][(*actHBC)->dof] = HDBC;
           bcCounter_[HDBC]++;
        }
      }

      //Get Grip of IdBC List for the fefunction
      const IdBcList idbcs = feFunction_->GetInHomDirichletBCs();
      IdBcList::const_iterator actIBC;

      for(actIBC = idbcs.Begin(); actIBC != idbcs.End(); actIBC++) {
        // Get all (Virtual) Nodes of the list
        GetNodesOfEntities(actNodes,(*actIBC)->entities);
        for(UInt iNode = 0 ; iNode < actNodes.GetSize();iNode++){
           //TODO find the source
           //make it zero based
           if( nodeMap_.BcKeys.find(actNodes[iNode]) == nodeMap_.BcKeys.end()){
             nodeMap_.BcKeys[actNodes[iNode]] = StdVector<BC_Type>(dofsPerUnknown,NOBC);
           }
           nodeMap_.BcKeys[actNodes[iNode]][(*actIBC)->dof] = IDBC;
           bcCounter_[IDBC]++;
        }
      }

      //Get Grip of constraint List for the fefunction
      const ConstraintList constraints = feFunction_->GetConstraints();
      ConstraintList::const_iterator actConstr;
      for(actConstr = constraints.Begin(); actConstr != constraints.End(); actConstr++) {
        StdVector<UInt> slaveNodes;
        GetNodesOfEntities(slaveNodes,(*actConstr)->slaveEntities);
        UInt masterDof = (*actConstr)->masterDof;
        UInt slaveDof = (*actConstr)->slaveDof;
        UInt mNode = slaveNodes[0];

        for ( UInt iNode = 1; iNode < slaveNodes.GetSize(); iNode++ ) {
          if( nodeMap_.BcKeys.find(slaveNodes[iNode]) == nodeMap_.BcKeys.end()){
            nodeMap_.BcKeys[slaveNodes[iNode]] = StdVector<BC_Type>(dofsPerUnknown,NOBC);
          }
          nodeMap_.BcKeys[slaveNodes[iNode]][slaveDof] = CONSTRAINT;
          nodeMap_.constraintNodes[std::pair<Integer,Integer>(slaveNodes[iNode],slaveDof)] = std::pair<Integer,Integer>(mNode,masterDof);
          bcCounter_[CONSTRAINT]++;
        }
      }

      //DEBUG output reaenable along with logging
      //std::map< Integer,StdVector<FeSpace::BC_Type> >::iterator iter = nodeMap_.BcKeys.begin();
      //for(iter; iter!= nodeMap_.BcKeys.end();iter++){
      //  std::cout << "The node #" << iter->first << " has the following flags: " << std::endl;
      //  for(UInt i = 0; i < iter->second.GetSize() ; i++){
      //    std::cout <<  iter->second[i] << std::endl;
      //  }
      //  std::cout << std::endl;
      //}
    }

    void FeSpaceH1Lagrange::CreateVirtualNodes(){
      //follow the following algorithm
      // - loop over every element
      //  - get edges and if not already mapped assign new virtual nodes
      //  - get faces and if not already mapped assign new virtual nodes
      //  - assign interior node to the element
      //  - now fill the Virtual node map with the information according to element orientation
      // - finally delete all intermediate arrays
      StdVector< shared_ptr<EntityList> > fctEntList;
      std::map<UInt, StdVector<Integer> > edgenodes;
      std::map<UInt, StdVector<Integer> > facenodes;
      std::map<UInt, StdVector<Integer> > interiornodes;

      fctEntList = feFunction_->GetEntityList();
      //get the highes possible node number
      UInt offset = feFunction_->GetGrid()->GetNumNodes();

      for(UInt actList = 0;actList <  fctEntList.GetSize(); actList++){
        // Get iterator of current element list
        EntityIterator entIt = fctEntList[actList]->GetIterator();


        //a little helper to create the nodes_ array
        if( mapType_ == GRID ) {
          StdVector<UInt> curNodes;
          GetNodesOfEntities( curNodes, fctEntList[actList] );
          for ( UInt aNode= 0; aNode < curNodes.GetSize(); aNode++ ) {
            if(gridToVirtualNodes_.find(curNodes[aNode]) == gridToVirtualNodes_.end()){
              gridToVirtualNodes_[curNodes[aNode]] =  curNodes[aNode];
            }
          }
        }

        for(entIt.Begin(); !entIt.IsEnd();entIt++){

          //check if we got what we expected
          if ( !entIt.GetType() == EntityList::ELEM_LIST ||
               !entIt.GetType() == EntityList::SURF_ELEM_LIST)  {
            EXCEPTION("FeSpaceH1Lagrange::CreateVirtualNodes(): Calling the method with an unsupported EntityListType");
          }

          const Elem* actEl = entIt.GetElem();
          BaseFE* ptFe = GetFe( entIt ); 

          //distinguish between Grid or polinomial based mapping
          if( mapType_ == GRID ) {
            StdVector<UInt> elemNodes = actEl->connect;
            virtualNodes_[actEl->elemNum][VERTEX] = elemNodes;
          } else if (mapType_ == POLYNOMIAL) {
            //get Grip of the reference element and set its order
            //ask FeFunction for the number of the order of approximation
            UInt approxOrder = feFunction_->GetIsoOrder();
            ptFe->SetIsoOrder(approxOrder);
            //===========================================================
            //Assign the Vertex node numbers
            //===========================================================
            StdVector<UInt> elemNodes = actEl->connect;
            UInt numVert = Elem::shapes[actEl->type].numVertices;
            virtualNodes_[actEl->elemNum][VERTEX].Resize(numVert);
            for ( UInt aNode= 0; aNode < numVert; aNode++ ) {
              virtualNodes_[actEl->elemNum][VERTEX][aNode] = elemNodes[aNode];
              if(gridToVirtualNodes_.find(elemNodes[aNode]) == gridToVirtualNodes_.end()){
                gridToVirtualNodes_[elemNodes[aNode]] =  elemNodes[aNode];
                nodes_.Push_back(elemNodes[aNode]);
              }
            }
            //Create the permutation array
            StdVector<UInt> permutations;

            if(approxOrder > 1){
              feFunction_->GetGrid()->MapEdges();
              feFunction_->GetGrid()->MapFaces();
            }
            //===========================================================
            //Assign the Edge node numbers
            //===========================================================
            UInt numEdgeNodes = 0;
            ElemShape actShape = Elem::shapes[actEl->type];
            for ( UInt iEdge=0; iEdge < actShape.numEdges; iEdge++) {
              UInt edgeNum = actEl->edges[iEdge];
              //get the permutation Vector
              ptFe->GetFncPermutation(permutations,actEl,EDGE,iEdge); 
              numEdgeNodes = permutations.GetSize(); 

              //Check if the current edge already has nodes or not
              if(edgenodes[edgeNum].GetSize() == 0){
                //here we assume spectral element approximation and we have 
                
                //order-1 nodes on the edge
                edgenodes[edgeNum].Resize(numEdgeNodes);
                for ( UInt edgeNode = 0;edgeNode < numEdgeNodes ;edgeNode++ ) {
                  edgenodes[edgeNum][edgeNode] = ++offset;
                  nodes_.Push_back(offset);
                }
              }

              //fill the virtual Nodes in the correct ordering
              for ( UInt i = 0; i < numEdgeNodes ; i++ ) {
                virtualNodes_[actEl->elemNum][EDGE].Push_back(edgenodes[edgeNum][ permutations[i] ]);
              }
            }

            //===========================================================
            //Assign the Face node numbers
            //===========================================================
            UInt numFaceNodes = 0; 
            for ( UInt iFace=0; iFace < actShape.numFaces; iFace++) {
              UInt faceNum = actEl->faces[iFace];
              //get the permutation Vector
              ptFe->GetFncPermutation(permutations,actEl,FACE,iFace); 
              numFaceNodes = permutations.GetSize(); 

              //Check if the current face already has nodes or not
              if(facenodes[faceNum].GetSize() == 0){
                //here we assume spectral element approximation and we have 
                //order-1 nodes on the edge
                facenodes[faceNum].Resize(numFaceNodes);
                for ( UInt faceNode = 0;faceNode < numFaceNodes ;faceNode++ ) {
                  facenodes[faceNum][faceNode] = ++offset;
                  nodes_.Push_back(offset);
                }
              }
              //fill the virtual Nodes in the correct ordering
              for ( UInt i = 0; i < numFaceNodes ; i++ ) {
                virtualNodes_[actEl->elemNum][FACE].Push_back(facenodes[faceNum][ permutations[i] ]);
              }
            }

            //===========================================================
            //Assign the Interior node numbers
            //===========================================================
            //get the permutation Vector just for the number of nodes
            ptFe->GetFncPermutation(permutations,actEl,INTERIOR,0); 
            UInt numIntNodes = permutations.GetSize(); 

            //Check if the current face already has nodes or not
            if(interiornodes[actEl->elemNum].GetSize() == 0){
              //here we assume spectral element approximation and we have 
              //order-1 nodes on the edge
              interiornodes[actEl->elemNum].Resize(numIntNodes);
              for ( UInt intNode = 0;intNode < numIntNodes ;intNode++ ) {
                facenodes[actEl->elemNum][intNode] = ++offset;
                nodes_.Push_back(offset);
              }
            }
            //fill the virtual Nodes in the correct ordering
            for ( UInt i = 0; i < numIntNodes ; i++ ) {
              virtualNodes_[actEl->elemNum][INTERIOR].Push_back(interiornodes[actEl->elemNum][ permutations[i] ]);
            }
          }
          
        }
      }
      //another little helper to create the nodes_ array
      if( mapType_ == GRID ) {
        nodes_.Resize(gridToVirtualNodes_.size());
        UInt counter = 0;
        for(std::map<UInt,UInt>::const_iterator it = gridToVirtualNodes_.begin(); it != gridToVirtualNodes_.end(); ++it) {
          nodes_[counter++] = it->second;
        }
      }
    }
}
