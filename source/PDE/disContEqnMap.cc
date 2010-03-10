#include "disContEqnMap.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "Domain/domain.hh" 
#include "Utils/coordSystem.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"

#include <boost/lexical_cast.hpp> 

using std::string; 
using boost::lexical_cast; 

namespace CoupledField {

  // declare class specific logging stream
  DECLARE_LOG(disContEqnMap)
  DEFINE_LOG(disContEqnMap, "disContEqnMap")


  DiscontinuousEqnMap::DiscontinuousEqnMap(Grid* grid, PdeIdType pdeId, bool usePenalty )//
                      :EqnMap(grid,pdeId,usePenalty) {

  }

  DiscontinuousEqnMap::DiscontinuousEqnMap(Grid* grid, PdeIdType pdeId, bool usePenalty,EqnMap* startMap) 
                      :EqnMap(grid,pdeId,usePenalty) {
    ptGrid_ = grid;

    usePenalty_ = usePenalty;
    pdeId_ = pdeId;
    isFinalized_ = false;


    numEqns_ = startMap->GetNumEqns();
    numIdBcs_ = startMap->GetNumInHomDirichletEqns();
    numCs_ = startMap->GetNumConstraintSlaveEqns();
    
    numLocNodes_ = 0;
    numLocElems_ = 0;
    numLocEdges_ = 0;
    numLocFaces_ = 0;
  }

  DiscontinuousEqnMap::~DiscontinuousEqnMap() {

  }

  //! ======================================================================
  //! SET AND INITIALIZATION METHODS
  //! ======================================================================

  void DiscontinuousEqnMap::ReorderMapping( const StdVector<UInt>& order ) {

    LOG_DBG(disContEqnMap) << "WARNING: The ReorderMapping functionality has never \
                         been tested for discontinuous equation numbering" 
                     << std::endl;
    
    // Iterate over all nodal maps and map the new ordering
    EqnMapType::iterator it1;

    // Check if any reordering array was given
    if( !order.GetSize() ) {
      //std::cerr << "Performing no reordering";
      return;
    }

    // iterate over all nodal resuls
    for ( it1=nodeEqns_.begin(); it1!=nodeEqns_.end(); it1++ ) {
      Matrix<Integer> & actMap = it1->second;

      // iterate over all entries in the current map and renumber them
      for ( UInt i = 0; i < actMap.GetNumRows(); i++ ) {
        for ( UInt j = 0; j < actMap.GetNumCols(); j++ ) {
          if ( actMap[i][j] > 0 ) {
            //std::cout << actMap[i][j]-1<< std::endl;
            actMap[i][j] = (Integer)order[actMap[i][j]-1];
          }
          else if(actMap[i][j] < 0 ) {
            //due to constraints
            actMap[i][j] = -(Integer)order[-actMap[i][j]-1];   
          }
        }
      }
    }
    
    VecEqnMapType::iterator it2;

    // iterate over all face
    for ( it2=faceEqns_.begin(); it2!=faceEqns_.end(); it2++ ) {
      StdVector<Vector <Integer> >& actMap = it2->second;
      
      // iterate over all entries in the current map and renumber them
      for ( UInt i = 0; i < actMap.GetSize(); i++ ) {
        for ( UInt j = 0; j < actMap[i].GetSize(); j++ ) {
          if ( actMap[i][j] > 0 ) {
            actMap[i][j] = order[actMap[i][j]-1];
          }
          else if(actMap[i][j] < 0 ) {
            //due to constraints
            actMap[i][j] = -order[-actMap[i][j]-1];   
          }
        }
      }
      
    }

    // iterate over all edge resuls
    for ( it2=edgeEqns_.begin(); it2!=edgeEqns_.end(); it2++ ) {
      StdVector<Vector <Integer> >& actMap = it2->second;
      
      // iterate over all entries in the current map and renumber them
      for ( UInt i = 0; i < actMap.GetSize(); i++ ) {
        for ( UInt j = 0; j < actMap[i].GetSize(); j++ ) {
          if ( actMap[i][j] > 0 ) {
            actMap[i][j] = order[actMap[i][j]-1];
          }
          else if(actMap[i][j] < 0 ) {
            //due to constraints
            actMap[i][j] = -order[-actMap[i][j]-1];   
          }
        }
      }
      
    }
    
    
    // Iterate over all element maps and map the new ordering
    for ( it2=elemEqns_.begin(); it2!=elemEqns_.end(); it2++ ) {
      StdVector<Vector < Integer> > & actMap = it2->second;
      
      // iterate over all entries in the current map and renumber them
      for ( UInt i = 0; i < actMap.GetSize(); i++ ) {
        for ( UInt j = 0; j < actMap[i].GetSize(); j++ ) {
          if ( actMap[i][j] > 0 ) {
            actMap[i][j] = order[actMap[i][j]-1];
          }
          else if(actMap[i][j] < 0 ) {
            //due to constraints
            actMap[i][j] = -order[-actMap[i][j]-1];   
          }
        }
      }
    }
  }


  
  void DiscontinuousEqnMap::GetEqns( StdVector<Integer>& eqns,
                        const ResultInfo& result, const EntityIterator& it,
                        UInt dof ) const 
  {

    // Note: this is currently hard-coded

    // first of all, delete eqns-array
    eqns.Clear();
    
    UInt numDofs = result.dofNames.GetSize();
    // ============
    //  NODAL PART 
    // ============

    if ( result.definedOn == ResultInfo::NODE ||
         result.definedOn == ResultInfo::PFEM  ) 
    {
      
      // get related nodal equaiton map
      Matrix<Integer> const & map = (nodeEqns_.find( result ) )->second;

        // Distinguish the type the of the list
        if ( it.GetType() == EntityList::ELEM_LIST ||
             it.GetType() == EntityList::SURF_ELEM_LIST ) 
        {
          StdVector<UInt> const & nodes = it.GetElem()->connect;

          eqns.Resize( nodes.GetSize());
          eqns.Init();
          for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) 
          {
            Integer localNode = -999;
              for (UInt i = 0 ;i < mesh2DisPdeNode_[ nodes[iNode]-1 ].GetSize() ; i++) {
                if(static_cast<UInt>(mesh2DisPdeNode_[ nodes[iNode]-1 ][i][0]) == it.GetElem()->elemNum)
                  localNode = mesh2DisPdeNode_[ nodes[iNode]-1 ][i][1];
              }

            if (localNode < 1 ) 
              eqns[iNode] = 0;
            else 
              eqns[iNode] = map[localNode-1][dof-1];
          }
        } else if ( it.GetType() == EntityList::NODE_LIST ) 
        {
          LOG_DBG3(disContEqnMap) << "WARNING: GetDisContEqn Called with nodelist. Its possible that this will return the wrong\
                      equation number! " << std::endl;
          UInt node = it.GetNode();
          eqns.Resize( 0);
          eqns.Init();
          
          Integer localNode = -999;
          eqns.Resize(0);
          for (UInt i = 0 ;i < mesh2DisPdeNode_[ node-1 ].GetSize() ; i++)
          {
            localNode = mesh2DisPdeNode_[ node-1 ][i][1];
            if (localNode < 1 ) {
              eqns.Push_back(0);
            } else {
              eqns.Push_back(map[localNode-1][dof-1]);
            }
          }
        } else 
        {
          EXCEPTION( "This type of entity list is not defined for " 
                       << "equation mapping!" );
        }
      }
    // ===============
    //   EDGE PART 
    // ===============
    
    if ( result.definedOn == ResultInfo::EDGE ||
        result.definedOn == ResultInfo::PFEM ) {
      
      // Check if entity is of type NODE:
      // In this case, we have to leave, as currently only nodal
      // boundary conditions can be assigned.
      if( it.GetType() == EntityList::NODE_LIST ) {
        return;
      }
      
      // get edge map of current result
      StdVector<Vector<Integer> > const & map = 
            (edgeEqns_.find( result ) )->second;
      
        // get edges
      StdVector<Integer> const & edges = it.GetElem()->edges;
      
      // iterate over all edges
      // This seems for now to be quite sensless commented out for testing purpose
      //if( it.GetType() == EntityList::SURF_ELEM_LIST && result.fctType->IsDiscontinuous()) 
      //{
      //  for( UInt iEdge = 0; iEdge < edges.GetSize(); iEdge++ ) 
      //  {
      //    Integer locEdge = -999;
      //    for (UInt i = 0 ;i < mesh2DisPdeEdge_[ std::abs(edges[iEdge])-1 ].GetSize() ; i++)
      //    {
      //      locEdge = mesh2DisPdeEdge_[ std::abs(edges[iEdge])-1 ][i][1];
      //      for( UInt iDof = dof-1 ; iDof < map[locEdge-1].GetSize(); iDof+=numDofs ) 
      //        eqns.Push_back( map[locEdge-1][iDof] );
      //    }
      //  }
      //}else
      for( UInt iEdge = 0; iEdge < edges.GetSize(); iEdge++ ) 
      {
        Integer locEdge = -999;
        for (UInt i = 0 ;i < mesh2DisPdeEdge_[ std::abs(edges[iEdge])-1 ].GetSize() ; i++)
        {
            if(static_cast<UInt>(mesh2DisPdeEdge_[ std::abs(edges[iEdge])-1 ][i][0]) == it.GetElem()->elemNum)
              locEdge = mesh2DisPdeEdge_[ std::abs(edges[iEdge])-1 ][i][1];
        }
  
        for( UInt iDof = dof-1 ; iDof < map[locEdge-1].GetSize(); iDof+=numDofs ){
          eqns.Push_back( map[locEdge-1][iDof] );
        }
      }
    }

    // ===============
    //   FACE PART 
    // ===============

    if ( result.definedOn == ResultInfo::FACE ||
          result.definedOn == ResultInfo::PFEM  ) 
    {
      // Check if entity is of type NODE:
      // In this case, we have to leave, as currently only nodal
      // boundary conditions can be assigned.
      if( it.GetType() == EntityList::NODE_LIST ) {
        return;
      }
      
      // get face map of current result
      StdVector<Vector<Integer> > const & map = (faceEqns_.find( result ) )->second;
      
      // get faces
      StdVector<Integer> const & faces = it.GetElem()->faces;
      
      // iterate over all faces
      for( UInt iFace = 0; iFace < faces.GetSize(); iFace++ ) 
      {
        // get local face number
        Integer locFace = 0;
        for (UInt i = 0 ;i < mesh2DisPdeFace_[ std::abs(faces[iFace])-1 ].GetSize() ; i++)
        {
            if(static_cast<UInt>(mesh2DisPdeFace_[ std::abs(faces[iFace])-1 ][i][0]) == it.GetElem()->elemNum)
              locFace = mesh2DisPdeFace_[ std::abs(faces[iFace])-1 ][i][1];
        }

        for( UInt iDof = dof-1 ; iDof < map[locFace-1].GetSize(); iDof+=numDofs ){ 
          eqns.Push_back( map[locFace-1][iDof] );
        }
      }
    }
        
    // =============
    //  ELEM PART
    // =============
    if( result.definedOn == ResultInfo::ELEMENT  ||
        result.definedOn == ResultInfo::PFEM ) 
    {
      StdVector< Vector<Integer> >const & elemMap = (elemEqns_.find( result ) )->second;

      Integer localElem = mesh2PdeElem_[(it.GetElem()->elemNum)-1];
      if (localElem < 1 ) {
         //  nothin to do here 
      } else 
      {
        eqns.Push_back( elemMap[localElem-1][dof-1] );
        LOG_DBG3(disContEqnMap) << "Pushin back eqn " <<  elemMap[localElem-1][dof-1]
                         << " for interior of element #" 
                         << it.GetElem()->elemNum << std::endl;
      }
    }

    LOG_DBG3(disContEqnMap) << "Equations are: " << eqns.Serialize();
    LOG_DBG3(disContEqnMap) << "Number of equations: " << eqns.GetSize() << std::endl;
  }
  
  void DiscontinuousEqnMap::GetEqns( StdVector<Integer> &eqns,
      const ResultInfo& result, 
      const EntityIterator& it) const{

    UInt numDofs = result.dofNames.GetSize();          

    // first of all, delete eqns-array
    eqns.Clear();
    // ============
    //  NODAL PART 
    // ============

    if ( result.definedOn == ResultInfo::NODE ||
          result.definedOn == ResultInfo::PFEM) {
      
      // get related nodal equaiton map
      Matrix<Integer> const & map = (nodeEqns_.find( result ) )->second;

      if(map.GetNumRows() !=0 && map.GetNumRows()!=0)
      {  
        // Distinguish the type the of the list
        if ( it.GetType() == EntityList::ELEM_LIST ||
              it.GetType() == EntityList::SURF_ELEM_LIST ) 
        {
          StdVector<UInt> const & nodes = it.GetElem()->connect;
          eqns.Resize( nodes.GetSize() * numDofs );
          eqns.Init();
  
          for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
            Integer localNode = -999;
            for (UInt i = 0 ;i < mesh2DisPdeNode_[ nodes[iNode]-1 ].GetSize() ; i++)
            {
              if(static_cast<UInt>(mesh2DisPdeNode_[ nodes[iNode]-1 ][i][0]) == it.GetElem()->elemNum)
                localNode = mesh2DisPdeNode_[ nodes[iNode]-1 ][i][1];
            }
    
            for (UInt iDof = 0; iDof < numDofs; iDof++ ) {
      
              if (localNode < 1 ) {
                eqns[iNode*numDofs + iDof] = 0;
              } else {
                eqns[iNode*numDofs + iDof] = map[localNode-1][iDof];
              }
            }
          }
        } else if ( it.GetType() == EntityList::NODE_LIST ) {
          UInt node = it.GetNode();
          UInt numDofs = result.dofNames.GetSize();
  
          Integer localNode = -999;
            eqns.Resize( numDofs*mesh2DisPdeNode_[ node-1 ].GetSize() );
            eqns.Init();
            for (UInt i = 0 ;i < mesh2DisPdeNode_[ node-1 ].GetSize() ; i++)
            {
              localNode = mesh2DisPdeNode_[ node-1 ][i][1];
              for (UInt iDof = 0; iDof < numDofs; iDof++ ) {
                if (localNode < 1 ) {
                  eqns[i*numDofs + iDof] = 0;
                } else {
                  eqns[i*numDofs + iDof] = map[localNode-1][iDof];
                }
              }
            }
        } else {
          EXCEPTION( "This type of entity list is not defined for " 
                   << "equation mapping!" );
        }
      }
    }

    // ===============
    //   EDGE PART 
    // ===============
    
    if ( result.definedOn == ResultInfo::EDGE ||
         result.definedOn == ResultInfo::PFEM) {
            
      // Check if entity is of type NODE:
      // In this case, we have to leave, as currently only nodal
      // boundary conditions can be assigned.
      if( it.GetType() == EntityList::NODE_LIST ) {
        return;
      }
      // get edge map of current result
      StdVector<Vector<Integer> > const & map = (edgeEqns_.find( result ) )->second;
      
      // get edges
      StdVector<Integer> const & edges = it.GetElem()->edges;
          
      // iterate over all edges
      for( UInt iEdge = 0; iEdge < edges.GetSize(); iEdge++ ) {
      
        // get local edge number
        Integer locEdge = -999;
        for (UInt i = 0 ;i < mesh2DisPdeEdge_[ std::abs(edges[iEdge])-1 ].GetSize() ; i++)
        {
          if(static_cast<UInt>(mesh2DisPdeEdge_[ std::abs(edges[iEdge])-1 ][i][0]) == it.GetElem()->elemNum)
            locEdge = mesh2DisPdeEdge_[ std::abs(edges[iEdge])-1 ][i][1];
        }
        
        // iterate over all dofs of this edge
        for( UInt iDof = 0; iDof < map[locEdge-1].GetSize(); iDof++ ) {
          eqns.Push_back( map[locEdge-1][iDof] );
        }
      }
    }  

    // ===============
    //   FACE PART 
    // ===============

    if ( result.definedOn == ResultInfo::FACE ||
          result.definedOn == ResultInfo::PFEM ) {
      
      // Check if entity is of type NODE:
      // In this case, we have to leave, as currently only nodal
      // boundary conditions can be assigned.
      if( it.GetType() == EntityList::NODE_LIST ) {
        return;
      }
      
      // get face map of current result
      StdVector<Vector<Integer> > const & map = (faceEqns_.find( result ) )->second;
      
      // get faces
      StdVector<Integer> const & faces = it.GetElem()->faces;
      
      // iterate over all faces
      for( UInt iFace = 0; iFace < faces.GetSize(); iFace++ ) {
  
        // get local face number
        Integer locFace = -999;
        for (UInt i = 0 ;i < mesh2DisPdeFace_[ std::abs(faces[iFace])-1 ].GetSize() ; i++)
        {
          if(static_cast<UInt>(mesh2DisPdeFace_[ std::abs(faces[iFace])-1 ][i][0]) == it.GetElem()->elemNum)
            locFace = mesh2DisPdeFace_[ std::abs(faces[iFace])-1 ][i][1];
        }
    
        // iterate over all dofs of this face
        for( UInt iDof = 0; iDof < map[locFace-1].GetSize(); iDof++ ) {
          eqns.Push_back( map[locFace-1][iDof] );
        }
      }
    }

    // =============
    //  ELEM PART
    // =============
    if( result.definedOn == ResultInfo::ELEMENT  ||
        result.definedOn == ResultInfo::PFEM ) 
    {
      StdVector< Vector<Integer> >const & elemMap = (elemEqns_.find( result ) )->second;
      
      Integer localElem = mesh2PdeElem_[(it.GetElem()->elemNum)-1];
      for (UInt iDof = 0; iDof < elemMap[localElem-1].GetSize(); iDof++ ) {
        if (localElem < 1 ) 
        {
          eqns.Push_back(0);
        } else 
        {
          eqns.Push_back( elemMap[localElem-1][iDof] );
          LOG_DBG3(disContEqnMap) << "Pushin back eqn " <<  elemMap[localElem-1][iDof]
                           << " for interior of element #" 
                           << it.GetElem()->elemNum << std::endl;
        }
      }
    }

    LOG_DBG3(disContEqnMap) << "Equations are: " << eqns.Serialize();
    LOG_DBG3(disContEqnMap) << "Number of equations: " << eqns.GetSize() << std::endl;
  }
    

  Integer DiscontinuousEqnMap::GetEqn( const ResultInfo& result, 
              const EntityIterator& it,
              UInt dof ) const
  {
    StdVector<Integer> eqns;
    GetEqns( eqns, result, it );
    return eqns[dof-1];
  }

    


  Integer DiscontinuousEqnMap::GetNodeEqn( const ResultInfo& result,
      UInt nodeNr, UInt dof ) { 
    
    EXCEPTION( "EXEPTION: This call is not well suited for discontinuous elements.\
                 USE:\n StdVector<Integer>  DiscontinuousEqnMap::GetNodeEqns( \
                 const ResultInfo& result, UInt nodeNr, UInt dof ) \n \
                 OR:\n \
                 Integer DiscontinuousEqnMap::GetNodeEqn( const ResultInfo& result,\
                 const ResultInfo& result, UInt ElemNum, UInt nodeNr, UInt dof )" );
    //if ( result.definedOn == ResultInfo::NODE
    //    || result.definedOn == ResultInfo::PFEM ) {
    //  
    //  Matrix<Integer> const & map = (nodeEqns_.find( result ) )->second;
    //  
    //  Integer localNode = mesh2PdeNode_[ nodeNr-1 ];
    //  
    //  if (localNode < 1 ) {
    //    return 0; 
    //  } else {
    //    return map[localNode-1][dof-1];
    //  }
    //  
    //}    
    //return 0;
  }

  Integer DiscontinuousEqnMap::GetNodeEqn( UInt nodeNr, UInt dof ) {

     EXCEPTION( "EXEPTION: This call is not well suited for discontinuous elements.\
                  USE:\n StdVector<Integer>  DiscontinuousEqnMap::GetNodeEqns( \
                  UInt nodeNr, UInt dof ) \n \
                  OR:\n \
                  Integer DiscontinuousEqnMap::GetNodeEqn( const ResultInfo& result,\
                  UInt ElemNum, UInt nodeNr, UInt dof )" );
    //// First check, if more than one type of results are defined
    //if ( nodeEqns_.size() != 1 ) {
    //  EXCEPTION( "GetNodeEqn() can only be used if exactly "
    //             << "one nodal result is present!" );
    //}

    //Matrix<Integer> const & map = (nodeEqns_.begin() )->second;
    //Integer localNode = mesh2PdeNode_[ nodeNr-1 ];
    //
    //if (localNode < 1 ) {
    //  return 0; 
    //} else {
    //  return map[localNode-1][dof-1];
    //}
  }

  
  void DiscontinuousEqnMap::GetNodeEqn( const StdVector<UInt>& nodeNrs, 
      StdVector<Integer>& eqnNrs ) {
     EXCEPTION( "EXEPTION: This call is not well suited for discontinuous elements.\
                  USE:\n void  DiscontinuousEqnMap::GetNodeEqns( const StdVector<UInt>& nodeNrs,\
                  StdVector< StdVector<Integer> >& eqnNrs ) \n \
                  OR:\n \
                  void DiscontinuousEqnMap::GetNodeEqn( UInt elemNum, const StdVector<UInt>& nodeNrs,\
                  StdVector< StdVector<Integer> >& eqnNrs )" );
    //// First check, if more than one type of results are defined
    //if ( nodeEqns_.size() != 1 ) {
    //  EXCEPTION( "GetNodeEqn() can only be used if exactly "
    //             << "one nodal result is present!" );
    //}
    //Matrix<Integer> const & map = (nodeEqns_.begin() )->second;
    //UInt dofsPerNode =  (nodeEqns_.begin() )->first.dofNames.GetSize();

    //eqnNrs.Resize( nodeNrs.GetSize() * dofsPerNode );
    //eqnNrs.Init();
    //
    //for( UInt iNode=0; iNode<nodeNrs.GetSize(); iNode++) {
    //  Integer localNode = mesh2PdeNode_[ nodeNrs[iNode]-1 ];
    //  
    //  for (UInt iDof = 0; iDof < dofsPerNode; iDof++ ) {

    //    if (localNode < 1 ) {
    //      eqnNrs[iNode*dofsPerNode + iDof] = 0;
    //    } else {
    //      eqnNrs[iNode*dofsPerNode + iDof] = map[localNode-1][iDof];
    //    }
    //  }
    //}
      
  }
   

    
  
  
  //! ======================================================================
  //! LOCAL/GLOBAL MAPPING OF MESH ENTITIES
  //! ======================================================================
  
  void DiscontinuousEqnMap::Mesh2PdeNode(StdVector<UInt> & PdeNodes,
      const StdVector<UInt> & MeshNodes) const {
    LOG_DBG(disContEqnMap) << "Mesh2PdeNode: Warining, just returns all equation numbers for the PDE without any special ordering"
             << std::endl;

    PdeNodes.Resize(numLocNodes_);
    PdeNodes.Init();
   
    for (UInt i=0; i<MeshNodes.GetSize(); i++){ 
      for (UInt i = 0 ;i < mesh2DisPdeNode_[ MeshNodes[i]-1 ].GetSize() ; i++) {
        PdeNodes[i] = mesh2DisPdeNode_[ MeshNodes[i]-1 ][i][1];
      }
    }
  }
 
  UInt DiscontinuousEqnMap::Mesh2PdeNode(const UInt meshNode) const {

     EXCEPTION( "EXEPTION: This call is not well suited for discontinuous elements.\
                  USE:\n StdVector<UInt>  DiscontinuousEqnMap::Mesh2PdeNode( const UInt meshNode ) \n \
                  OR:\n \
                  UInt DiscontinuousEqnMap::Mesh2PdeNode( UInt elemNum, const UInt meshNode )" );
    //if ( mesh2PdeNode_[meshNode-1] < 0 ) {
    //  EXCEPTION(  "MeshNode Nr. " << meshNode 
    //              << " has no local node number!" );
    //}

    //return abs(mesh2PdeNode_[meshNode-1]);
  }
    
  //! ======================================================================
  //! MISCELLANEOUS
  //! ======================================================================
  
  void DiscontinuousEqnMap::ToInfo(InfoNode* base) const
  {
    InfoNode* lg = base->Get("DiscontinuousLocalGlobal");

    // local <-> global mapping (nodes)
    InfoNode* root = lg->Get("nodes");
    root->Get("globalNodes")->SetValue(ptGrid_->GetNumNodes());
    root->Get("localNodes")->SetValue(numLocNodes_);
    for( UInt iNode = 0; iNode < pde2MeshNode_.GetSize(); iNode++ )
    {
      InfoNode* in = root->Get("mapping", InfoNode::APPEND);
      in->Get("local")->SetValue(iNode+1);
      in->Get("global")->SetValue(pde2MeshNode_[iNode]);
    }

    // local <-> global mapping (elements)
    root = lg->Get("elements");
    root->Get("globalElements")->SetValue(ptGrid_->GetNumVolElems());
    root->Get("localElements")->SetValue(numLocElems_);
    for( UInt iElem = 0; iElem < pde2MeshElem_.GetSize(); iElem++ )
    {
      InfoNode* in = root->Get("mapping", InfoNode::APPEND);
      in->Get("local")->SetValue(iElem+1);
      in->Get("global")->SetValue(pde2MeshElem_[iElem]);
    }

    // local <-> global mapping (edges)
    root = lg->Get("edges");
    for( UInt iEdge = 0; iEdge < mesh2DisPdeEdge_.GetSize(); iEdge++ )
    {
      for( UInt edgeIdx = 0; edgeIdx < mesh2DisPdeEdge_[iEdge].GetSize(); edgeIdx++)
      {
        InfoNode* in = root->Get("mapping", InfoNode::APPEND);
        in->Get("local")->SetValue(iEdge+1);
        in->Get("global")->SetValue(mesh2DisPdeEdge_[iEdge][edgeIdx][1]);
      }
    }

    // local <-> global mapping (faces)
    root = lg->Get("faces");
    for( UInt iFace = 0; iFace < mesh2DisPdeFace_.GetSize(); iFace++ )
    {
      for( UInt faceIdx = 0;faceIdx<  mesh2DisPdeFace_[iFace].GetSize(); faceIdx++)
      {
        InfoNode* in = root->Get("mapping", InfoNode::APPEND);
        in->Get("local")->SetValue(iFace+1);
        in->Get("global")->SetValue(mesh2DisPdeFace_[iFace][faceIdx][1]);
      }
    }

    // Nodal equation mapping
    InfoNode* em = base->Get("DiscontinuousEquationMapping");
    InfoNode* rt = em->Get("nodal");

    // Loop over all nodal mapped results
    ResultEntityMap::const_iterator listIt;
    for(listIt = nodeMappedList_.begin(); listIt != nodeMappedList_.end();  listIt++)
    {
      // get dofspernode and associated eqnMap
      Matrix<Integer> const & eqnMap = (nodeEqns_.find(listIt->first))->second;
      UInt dofsPerNode = listIt->first.dofNames.GetSize();

      // head for which result this is shown
      InfoNode* res = rt->Get("result", InfoNode::APPEND);
      res->Get("type")->SetValue(SolutionTypeEnum.ToString(listIt->first.resultType));
      res->Get("dofs")->SetValue(dofsPerNode);

      for ( UInt iNode = 0; iNode < pde2MeshNode_.GetSize(); iNode++ )
      {
        InfoNode* in = res->Get("mapping", InfoNode::APPEND);
        in->Get("local")->SetValue(iNode+1);
        in->Get("global")->SetValue(pde2MeshNode_[iNode]);

        for ( UInt iDof = 0; iDof < dofsPerNode; iDof++ )
          in->Get("eqnNr_dof" + lexical_cast<string>(iDof+1))->SetValue(eqnMap[iNode][iDof]);
      }
    }

    // Loop over all edge mapped results
    rt = em->Get("edge");

    ResultEntityMap::const_iterator edgeListIt;
    for(edgeListIt = edgeMappedList_.begin(); edgeListIt != edgeMappedList_.end(); edgeListIt++ )
    {
      StdVector<Vector<Integer> >  const & eqnMap =
        (edgeEqns_.find(edgeListIt->first))->second;

      // Print head for which result this is shown
      InfoNode* res = rt->Get("result", InfoNode::APPEND);
      res->Get("type")->SetValue(SolutionTypeEnum.ToString(edgeListIt->first.resultType));

      for ( UInt iEdge = 0; iEdge < ptGrid_->GetNumEdges(); iEdge++ )
      {
        for( UInt edgeIdx = 0; edgeIdx < mesh2DisPdeEdge_[iEdge].GetSize(); edgeIdx++)
        {
          // local edge numger
          Integer locEdge = mesh2DisPdeEdge_[iEdge][edgeIdx][1];
          // check if edge has any equations at all
          if ( locEdge <= 0 ) break;

          InfoNode* in = res->Get("mapping", InfoNode::APPEND);
          in->Get("local")->SetValue(iEdge+1);
          in->Get("global")->SetValue(mesh2DisPdeEdge_[iEdge][edgeIdx][1]);

          if(eqnMap[locEdge-1].GetSize() > 0)
          {
            for (UInt iDof = 0; iDof < eqnMap[locEdge-1].GetSize(); iDof++ )
              in->Get("eqnNr_dof" + lexical_cast<string>(iDof + 1))->SetValue(eqnMap[locEdge-1][iDof]);
          }
        }
      }
    }

    // --- Face equation mapping ---
    rt = em->Get("face");

    // Loop over all face mapped results
    ResultEntityMap::const_iterator faceListIt;
    for(faceListIt = faceMappedList_.begin(); faceListIt != faceMappedList_.end(); faceListIt++)
    {
      StdVector<Vector<Integer> >  const & eqnMap =
        (faceEqns_.find(faceListIt->first))->second;

      // Print head for which result this is shown
      InfoNode* res = rt->Get("result", InfoNode::APPEND);
      res->Get("type")->SetValue(SolutionTypeEnum.ToString(faceListIt->first.resultType));

      for ( UInt iFace = 0; iFace < ptGrid_->GetNumFaces(); iFace++ )
      {
        for( UInt faceIdx = 0; faceIdx < mesh2DisPdeFace_[iFace].GetSize(); faceIdx++)
        {
          // local face numer
          Integer locFace = mesh2DisPdeFace_[iFace][faceIdx][1];

          // check if face has any equations at all
          if ( locFace <= 0 ) break;

          InfoNode* in = res->Get("mapping", InfoNode::APPEND);
          in->Get("localFaceNr")->SetValue(iFace+1);
          in->Get("globalFaceNr")->SetValue(mesh2DisPdeFace_[iFace][faceIdx][1]);

          if(eqnMap[locFace-1].GetSize() > 0)
          {
            for (UInt iDof = 0; iDof < eqnMap[iFace].GetSize(); iDof++ )
              in->Get("eqnNr_dof" + lexical_cast<string>(iDof+1))->SetValue(eqnMap[locFace-1][iDof]);
          }
        }
      }
    }

    // --- Interior equation mapping ---
    rt = em->Get("interior");

    // Loop over all element interior mapped results
    ResultEntityMap::const_iterator it;
    for(it = elemIntMappedList_.begin(); it != elemIntMappedList_.end(); ++it)
    {
      StdVector<Vector<Integer> >  const & eqnMap =
        (elemEqns_.find(it->first))->second;

      // Print head for which result this is shown
      InfoNode* res = rt->Get("result", InfoNode::APPEND);
      res->Get("type")->SetValue(SolutionTypeEnum.ToString(it->first.resultType));

      for ( UInt iElem = 0; iElem < numLocElems_; iElem++ )
      {
        InfoNode* in = res->Get("mapping", InfoNode::APPEND);
        in->Get("localElemNr")->SetValue(iElem+1);
        in->Get("globalElenNr")->SetValue(pde2MeshElem_[iElem]);

        if(eqnMap[iElem].GetSize() > 0)
        {
          for (UInt iDof = 0; iDof < eqnMap[iElem].GetSize(); iDof++ )
            in->Get("eqnNr_dof" + lexical_cast<string>(iDof+1))->SetValue(eqnMap[iElem][iDof]);
        }
        else in->Get("eqnNr_dof1")->SetValue(0);
      }
    } 
  }


  //! ======================================================================
  //! PRIVATE HELPER METHODS
  //! ======================================================================

  
  void DiscontinuousEqnMap::CalcNodeElemMapping() {
    
    mesh2DisPdeNode_.Resize( ptGrid_->GetNumNodes() );
    pde2MeshNode_.Clear( );

    mesh2PdeElem_.Resize( ptGrid_->GetNumElems() );
    mesh2PdeElem_.Init( -1 );
    // Note: here we could iterate over all element lists and 
    //       count the number of entries, so we would now from beginning,
    //       how many elements we have in this pde
    pde2MeshElem_.Clear( );
 

    //UInt nodeCounter = 0;
    UInt elemCounter = 1;
    StdVector<Elem*> subdom;
 
    // Note: The check, if the elemList of this region corresponds to all 
    //       elements of the pde can not yet be performed
    // First, check, if Pde is defined on all regions of domain
    //  StdVector<RegionIdType> allRegions;
    bool allFound = false;
    //     Integer index = -1;

    //     ptGrid_->GetVolRegionIds( allRegions );
    //     for ( UInt i = 0; i < allRegions.GetSize(); i++ ) {
    //       index = subdoms_.Find(allRegions[i]);
    //       if ( index == -1 ) {
    //         allFound = false;
    //         break;
    //       }
    //     }
    
    // Case 1: This pde is defined on all volume regions of the grid
    //         --> Do a simple 1-to-1 mapping of global node number
    //             and local one
    if ( allFound ) {
      //    pde2MeshNode_.Resize(mesh2PdeNode_.GetSize());
      //       for (UInt i = 0; i<mesh2PdeNode_.GetSize(); i++) {
      //         mesh2PdeNode_[i] = i+1;
      //         pde2MeshNode_[i] = i+1;
      //     }
      
      //       pde2MeshElem_.Resize(ptGrid_->GetNumElems());
      //       for (UInt i = 0; i<mesh2PdeElem_.GetSize(); i++) {
      //         pde2MeshElem_[i] = i+1;
      //         mesh2PdeElem_[i] = i+1;
      //       }
      
      
      //       nodeCounter=mesh2PdeNode_.GetSize();
  } else 
  {
      // Case 2: This pde is defined on a subset of all regions
      //         --> Perform normal node renumbering
      
      
      //  // iterate over all element lists
    for ( UInt iList = 0; iList < locEntities_.GetSize(); iList++ ) 
    {
      // Get iterator of current element list
      EntityIterator it = locEntities_[iList]->GetIterator();
     
      // 1) Check, if entity list contains elements
      if( locEntities_[iList]->GetType() == EntityList::ELEM_LIST ||
          locEntities_[iList]->GetType() == EntityList::SURF_ELEM_LIST ) {
        // iterate over all elements in element list
        for ( it.Begin(); !it.IsEnd(); it++ ) {
          
          // Store current element
          const Elem* actEl = it.GetElem();
          // *** Mapping of Elements ***
          mesh2PdeElem_[actEl->elemNum - 1 ] = elemCounter;
          pde2MeshElem_.Push_back( actEl->elemNum );
          elemCounter++;
        }
      }
      
      // 2) Perform nodal mapping in any case
      it = locEntities_[iList]->GetIterator();
      for ( it.Begin(); !it.IsEnd(); it++ ) 
      {
        const Elem* actEl = it.GetElem();
        
        //this check should not be necessary happen anyway
        for ( UInt iNode=0; iNode < actEl->connect.GetSize(); iNode++)
        {
          bool found = false;
          for ( UInt i =0;i<mesh2DisPdeNode_[actEl->connect[iNode]-1].GetSize() ;i++ )
          {
            if(static_cast<UInt>(mesh2DisPdeNode_[actEl->connect[iNode]-1][i][0]) == actEl->elemNum)
            {  
              found = true;
              //std::cerr << "Error while mapping Node #"<< actEl->connect[iNode]-1 << "of element #" <<  actEl->elemNum << std::endl;
              //EXCEPTION("found an already mapped entity, this should not happen!");
              break;
            }
          }
          if(!found)
          {
            mesh2DisPdeNode_[actEl->connect[iNode]-1].Push_back(new Integer[2]);
            mesh2DisPdeNode_[actEl->connect[iNode]-1][mesh2DisPdeNode_[actEl->connect[iNode]-1].GetSize()-1][0] = actEl->elemNum;
            mesh2DisPdeNode_[actEl->connect[iNode]-1][mesh2DisPdeNode_[actEl->connect[iNode]-1].GetSize()-1][1] = ++numLocNodes_;
            pde2MeshNode_.Push_back(actEl->connect[iNode]);
            //std::cout << "Adding node " << actEl->connect[iNode] << " to element #" << mesh2DisPdeNode_[actEl->connect[iNode]-1][mesh2DisPdeNode_[actEl->connect[iNode]-1].GetSize()-1][0] << " with local node" << numLocNodes_ << std::endl;
          }
        }
      }
    }
  }

  // remember number of local nodes and elements
  numLocNodes_ = pde2MeshNode_.GetSize();
  numLocElems_ = pde2MeshElem_.GetSize();
}
  
  
  void DiscontinuousEqnMap::CalcEdgeMapping()  {

    mesh2DisPdeEdge_.Resize( ptGrid_->GetNumEdges() );
    
    // iterate over all element lists
    for ( UInt iList = 0; iList < locEntities_.GetSize(); iList++ ) {
      if( locEntities_[iList]->GetType() == EntityList::ELEM_LIST ||
          locEntities_[iList]->GetType() == EntityList::SURF_ELEM_LIST ) {
        
        EntityIterator it = locEntities_[iList]->GetIterator();
        // iterate over all elements in element list
        for ( it.Begin(); !it.IsEnd(); it++ ) 
        {
          // Store current element
          const Elem* actEl = it.GetElem();
          // iterate over all edges
          for ( UInt iEdge=0; iEdge < actEl->edges.GetSize(); iEdge++) 
          {
            bool found = false;
            for ( UInt i =0;i<mesh2DisPdeEdge_[std::abs(actEl->edges[iEdge])-1].GetSize() ;i++ )
            {
              if(static_cast<UInt>(mesh2DisPdeEdge_[std::abs(actEl->edges[iEdge])-1][i][0]) == actEl->elemNum)
              {  
                found = true;
                break;
              }
            }
            if(!found)
            {
              mesh2DisPdeEdge_[std::abs(actEl->edges[iEdge])-1].Push_back(new Integer[2]);
              mesh2DisPdeEdge_[std::abs(actEl->edges[iEdge])-1]
                              [mesh2DisPdeEdge_[std::abs(actEl->edges[iEdge])-1].GetSize()-1][0] 
                              = actEl->elemNum;
              mesh2DisPdeEdge_[std::abs(actEl->edges[iEdge])-1]
                              [mesh2DisPdeEdge_[std::abs(actEl->edges[iEdge])-1].GetSize()-1][1] 
                              = ++numLocEdges_;
            }
          }
        }
      }
    }
  }
  

  void DiscontinuousEqnMap::CalcFaceMapping()  {

    LOG_TRACE(disContEqnMap) << "Starting local<->global face mapping\n";

    mesh2DisPdeFace_.Resize( ptGrid_->GetNumFaces() );

    //UInt faceCounter = 0;

    // iterate over all element lists
    for ( UInt iList = 0; iList < locEntities_.GetSize(); iList++ ) {
      if( locEntities_[iList]->GetType() == EntityList::ELEM_LIST ||
          locEntities_[iList]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EntityIterator it = locEntities_[iList]->GetIterator();
        // iterate over all elements in element list
        for ( it.Begin(); !it.IsEnd(); it++ ) 
        {
          // Store current element
          const Elem* actEl = it.GetElem();
          // iterate over all faces
          for ( UInt iFace=0; iFace < actEl->faces.GetSize(); iFace++) 
          {
            bool found = false;
            for ( UInt i =0;i<mesh2DisPdeFace_[std::abs(actEl->faces[iFace])-1].GetSize() ;i++ )
            {
              if(static_cast<UInt>(mesh2DisPdeFace_[std::abs(actEl->faces[iFace])-1][i][0]) == actEl->elemNum)
              {  
                found = true;
                break;
              }
            }
            if(!found)
            {
              mesh2DisPdeFace_[std::abs(actEl->faces[iFace])-1].Push_back(new Integer[2]);
              mesh2DisPdeFace_[std::abs(actEl->faces[iFace])-1][mesh2DisPdeFace_[std::abs(actEl->faces[iFace])-1].GetSize()-1][0] = actEl->elemNum;
              mesh2DisPdeFace_[std::abs(actEl->faces[iFace])-1][mesh2DisPdeFace_[std::abs(actEl->faces[iFace])-1].GetSize()-1][1] = ++numLocFaces_;
            }//element not in list? 
          }//over all faces (for-loop)
        }//over all elements (for-loop)
      }//if theres an element
    }//over all element-lists (for-loop)

    LOG_DBG(disContEqnMap) << "There are " << numLocFaces_ << " local faces";
    LOG_TRACE(disContEqnMap) << "Finished local<->global face mapping\n";

  }
  
  void DiscontinuousEqnMap::CalcNodalEquations( UInt phase ) {

    // MAGIC number, which gets assignetd to all nodes,
    // which have not yet an equation number
    const Integer NO_EQN = -333;

    // Big outer loop over all nodal mapped element lists
    ResultEntityMap::iterator listIt;

    for( listIt = nodeMappedList_.begin(); 
    listIt != nodeMappedList_.end(); 
    listIt++ ) {

      // Remeber current result and list of elementLists
      const ResultInfo & actRes = listIt->first;
      StdVector<shared_ptr<EntityList> > & actLists = listIt->second;


      // Get grip of homogeneous and in-homogeneous boundary conditions
      // for this tpye of result
      ResultHdBcMap::iterator hdBcIt = hdBcs_.find( actRes );
      ResultIdBcMap::iterator idBcIt = idBcs_.find( actRes );
      ResultConstraintMap::iterator csIt = constraints_.find( actRes );

      Matrix<Integer> & actMap = nodeEqns_[actRes];

      // Get number of dofs
      UInt dofsPerNode = actRes.dofNames.GetSize();
    
      // Idea of the algorithm:
      //
      // -- PHASE 1 --
      // Step 1:  Initialize actMap with NO_EQN
      // Step 2:  For each entry in constraintSlaveNodes_ set the corresponding
      //          entry in pdeNode2eqn_ to 0
      // Step 3:  For each entry in homoDirichletNodes_ set the corresponding
      //          entry in pdeNode2eqn_ to 0
      // Step 3b: For each entry in inhomoDirichletNodes_ set the corresponding
      //          entry in pdeNode2eqn_ to -1
      
      // Step 4:  Loop over all nodes of given entity lists
      //          and assign each non-zero entry an equation number
      // Step 5:  Loop over the whole map and set all entries with NO_EQN to 0

      // -- PHASE 2 --
      // Step 6: Loop over all entries in inhomoDirichletNodes_ and assign that
      //          dof an equation number after the hightest equation number of
      //          the free dofs
      // Step 6b:  Afterwards loop again over all nodes in constraintSlaveNodes_
      //          and set the corresponding entry in pdeNode2EQN_ to the
      //          negative of the value of constraintMasterNode
      //
      
      if( phase == 1 ) {

        // ------
        // STEP 1
        // ------
        //UInt multipleBCs = 0;

        actMap.Resize( numLocNodes_ , dofsPerNode);
        actMap.InitValue( NO_EQN );

        // ------
        // STEP 2
        // ------

        // Check if any inhom. constraint is defined for the current
        // result
        if( csIt != constraints_.end() ) {
          ConstraintList const & actCsList = csIt->second;

          for ( UInt i = 0; i < actCsList.GetSize(); i++ ) {
            StdVector<UInt> slaveNodes;
            GetNodesOfEntities( slaveNodes, actCsList[i]->slaveEntities );
            UInt slaveDof = actCsList[i]->slaveDof;

            // NOTE: at the moment we assume that slave and master nodes
            // are the same, therefore we start with the second node
            // within the master / slave node array
            for ( UInt iNode = 1; iNode < slaveNodes.GetSize(); iNode++ ) {
                for ( UInt i=0;i<mesh2DisPdeNode_[slaveNodes[iNode]-1].GetSize() ;i++ )
                  actMap[mesh2DisPdeNode_[slaveNodes[iNode]-1][i][1]-1][slaveDof-1] = -1;
            }
          }
        }


        // ------
        // STEP 3
        // ------
        Matrix<UInt> countNodes;
        countNodes.Resize( numLocNodes_, dofsPerNode );
        countNodes.Init();

        if( hdBcIt != hdBcs_.end() ) {
          HdBcList const & actHdBcList = hdBcIt->second;
          for ( UInt i = 0; i < actHdBcList.GetSize(); i++ ) {
            StdVector<UInt> nodes;
            GetNodesOfEntities( nodes, actHdBcList[i]->entities );
            UInt actDof = actHdBcList[i]->dof;

            for( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
              // Check if homDirichletNode belongs to one of my subdomains
              if ( mesh2DisPdeNode_[nodes[iNode]-1].GetSize() == 0 ) {
                (*warning) << "DiscontinuousEqnMap::CalcNodalEquations: Homogen. Dirichlet node "
                << "nr. " << nodes[iNode]
                                   << " is not contained in any of the regions for this PDE";
                Warning( __FILE__, __LINE__ );
              }
              //else if ( countNodes[mesh2PdeNode_[nodes[iNode]-1]-1] [actDof-1] != 0 ) {
              //  //     (*warning) << "DiscontinuousEqnMap::CalcNodalEquations: HomDirichletNode # "
              //  //                << nodes[i]
              //  //                << "\nappeared already at least once in the list of "
              //  //                << "boundary nodes for this PDE!\n Please check, if this "
              //  //                << "node is defined in more than one level of boundary "
              //  //                << "nodes!";
              //  //     Warning( __FILE__, __LINE__ );
              //}
              else {
                for ( UInt i=0;i<mesh2DisPdeNode_[nodes[iNode]-1].GetSize() ;i++ )
                {
                  if ( countNodes[mesh2DisPdeNode_[nodes[iNode]-1][i][1]-1][actDof-1] != 0 )
                    continue;

                  actMap[mesh2DisPdeNode_[nodes[iNode]-1][i][1]-1][actDof-1] = 0;
                  countNodes[mesh2DisPdeNode_[nodes[iNode]-1][i][1]-1][actDof-1]++;
                }
              }
            }
          }
        }

        // -------
        // STEP 3b
        // -------

        countNodes.Init();

        // Check if any inhom. boundary condition is defined for the current
        // result
        if( idBcIt != idBcs_.end() ) {
          IdBcList const & actIdBcList = idBcIt->second;

          for ( UInt i = 0; i < actIdBcList.GetSize(); i++ ) {
            StdVector<UInt> nodes;
            GetNodesOfEntities( nodes, actIdBcList[i]->entities );
            UInt actDof = actIdBcList[i]->dof;

            for( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {

              if ( mesh2DisPdeNode_[ nodes[iNode] - 1 ].GetSize() == 0 ) {
                (*warning) << "DiscontinuousEqnMap::CalcNodalEquations: Inhom. Dirichlet "
                << "node #" << nodes[iNode]
                                     << " is not contained in any of the regions for "
                                     << "this Pde";
                Warning( __FILE__, __LINE__ );
              }
              //else if ( countNodes[mesh2PdeNode_[nodes[iNode]-1]-1]
              //                     [actDof-1] != 0 ) {
              //  //   (*warning) << "DiscontinuousEqnMap::CalcNodalEquations: Inhom. Dirichlet "
              //  //                << "node #" << nodes[iNode]
              //  //                << "\nappeared already at least once in the list of "
              //  //                << "boundary nodes for this Pde!\n Please check, if "
              //  //                << "this node is defined in more than one level of "
              //  //                << "boundary nodes!";
              //  //     Warning( __FILE__, __LINE__ );
              //}
              else {
                
                // only set entry to -1, if entry is not yet an constraint
                // slave entry or homogeneous dirichlet entry
                for ( UInt i=0;i<mesh2DisPdeNode_[nodes[iNode]-1].GetSize() ;i++ )
                {
                  if(  actMap[mesh2DisPdeNode_[nodes[iNode]-1][i][1]-1] [actDof-1] == NO_EQN ) {
                    if ( countNodes[mesh2DisPdeNode_[nodes[iNode]-1][i][1]-1][actDof-1] != 0 )
                      continue;
                    actMap[mesh2DisPdeNode_[nodes[iNode]-1][i][1]-1][actDof-1] = -1;
                    countNodes[mesh2DisPdeNode_[nodes[iNode]-1][i][1]-1][actDof-1]++;
                    // In any case we have to increment the number of idBC-conditions
                    numIdBcs_++;
                  }
                } 
              }
            }
          }
        }

        // ------
        // STEP 4
        // ------
        // Initialize countNodes to zero. It will be used to count if
        // a node got already an equation number
        countNodes.Init();

        StdVector<UInt> nodes;
        // Iterate over all element list belonging to this result
        // and assign each of them a node number
        for( UInt iList = 0; iList < actLists.GetSize(); iList++ ) {

          // Get nodes of current entityList

          Integer locNode = 0;
					//gehe über jedes element, hole knoten und nummeriere wenn elementNum in mesh2DisPdeNode_
					EntityIterator entIt = actLists[iList]->GetIterator();
					for(entIt.Begin(); !entIt.IsEnd(); entIt++)
					{
          	const Elem* actEl = entIt.GetElem();
						for(UInt iNode = 0; iNode < actEl->connect.GetSize(); iNode++)
						{
          		for ( UInt i=0;i<mesh2DisPdeNode_[actEl->connect[iNode]-1].GetSize() ;i++ )
          		{
								if(static_cast<UInt>(mesh2DisPdeNode_[actEl->connect[iNode]-1][i][0]) == actEl->elemNum) {
          		  	locNode = mesh2DisPdeNode_[actEl->connect[iNode]-1][i][1];
          		  	for ( UInt iDof = 0; iDof < dofsPerNode; iDof++ ) 
          		  	{
          		  	  if ( actMap[locNode-1][iDof] == NO_EQN )//&& countNodes[locNode-1][iDof] == 0 ) 
          		  	  {
          		  	    numEqns_++;
          		  	    LOG_DBG3(disContEqnMap) << "Adding equation " << numEqns_ 
          		  	                     << " to local node " << locNode << std::endl;
          		  	    //std::cout << "Adding equation " << numEqns_ 
          		  	    //                 << " to local node " << locNode << " to DOF" << iDof << std::endl;
          		  	    actMap[locNode-1][iDof] = numEqns_;
          		  	    countNodes[locNode-1][iDof] += 1;
          		  	  }
          		  	}
								}
          		}
						}
					}
        }

        // ------
        // STEP 5
        // ------
        // Re-iterate over the whole equation map and set all entries 
        // with eqn-number of NO_EQN to 0
        for( UInt i = 0; i < actMap.GetNumRows(); i++ ) {
          for( UInt j = 0; j < actMap.GetNumCols(); j++ ) {
            if( actMap[i][j] == NO_EQN )
              actMap[i][j] = 0;
          }
        }

        LOG_DBG2(disContEqnMap) << "Final equation map looks like: \n" 
                         << actMap << std::endl;
        
      } else if( phase == 2 ) {

        Integer locNode = 0;
        // -------
        // STEP 6
        // -------
        if( idBcIt != idBcs_.end() ) {
          IdBcList const & actIdBcList = idBcIt->second;
          for ( UInt i = 0; i < actIdBcList.GetSize(); i++ ) {
            StdVector<UInt> nodes;
            GetNodesOfEntities( nodes, actIdBcList[i]->entities );
            UInt actDof = actIdBcList[i]->dof;

            for ( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
              for ( UInt i=0;i<mesh2DisPdeNode_[nodes[iNode]-1].GetSize() ;i++ )
              {
                locNode = mesh2DisPdeNode_[nodes[iNode]-1][i][1];
                if( locNode > 0 ) 
                { 
                  numEqns_++;
                  actMap[locNode-1] [actDof-1] = numEqns_;
                }
              }
            }
          }
        }

        // ------
        // STEP 6b
        // ------
        // Check if any inhom. constraint is defined for the current
        // result
        if( csIt != constraints_.end() ) {
          ConstraintList const & actCsList = csIt->second;

          for ( UInt i = 0; i < actCsList.GetSize(); i++ ) {
            StdVector<UInt> slaveNodes;
            GetNodesOfEntities( slaveNodes, actCsList[i]->slaveEntities );
            UInt slaveDof = actCsList[i]->slaveDof;
            UInt masterDof = actCsList[i]->masterDof;
            UInt masterNode = slaveNodes[0];

            for ( UInt iNode = 1; iNode < slaveNodes.GetSize(); iNode++ ) 
            {
              for ( UInt i = 0;i<mesh2DisPdeNode_[slaveNodes[iNode]-1].GetSize() ;i++ )
              {
                actMap[mesh2DisPdeNode_[slaveNodes[iNode]-1][i][1]-1] [slaveDof-1] = 
                    -actMap[mesh2DisPdeNode_[masterNode-1][i][1]-1] [masterDof-1];

                numCs_++;
              }
            }
          }
        }

      } else {
        EXCEPTION( "Phase '" << phase << "' does not exist!" );
      }
    }
  }
  
  void DiscontinuousEqnMap::CalcElemInteriorEquations( UInt phase ) {
    
    // Big outer loop over all element interior mapped lists
    ResultEntityMap::iterator listIt;
    
    for( listIt = elemIntMappedList_.begin(); 
         listIt != elemIntMappedList_.end(); 
         listIt++ ) {

      // Remeber current result and list of elementLists
      const ResultInfo & actRes = listIt->first;
      StdVector<shared_ptr<EntityList> > & actLists = listIt->second;
      StdVector<Vector<Integer> > & actMap = elemEqns_[actRes];
      
      UInt dofsPerElem = actRes.dofNames.GetSize();
      actMap.Resize( numLocElems_  );

      // Fetch also nodal map
      Matrix<Integer> & actNodeMap = nodeEqns_[actRes];

      Integer locElem = 0;

      // Iterate over all entitylists
      for( UInt iList = 0; iList < actLists.GetSize(); iList++ ) {

        EntityIterator it = actLists[iList]->GetIterator();

        // Iterate over all elements within this list
        for( it.Begin(); !it.IsEnd(); it++ ) {

          // Get grip of element
          const Elem & actEl = *(it.GetElem());

          // Get number of unknowns for element
          StdVector< StdVector<UInt> > numFcns; 
          numFcns.Resize( dofsPerElem );

          // iterate over all dofs of this result
          for( UInt iDof = 0; iDof < dofsPerElem; iDof++ ) {
            actEl.ptElem->GetNumFncs( numFcns[iDof], actRes.fctType, 
                                      AnsatzFct::INTERIOR, iDof );
          }

          // determine local element number
          locElem = mesh2PdeElem_[actEl.elemNum - 1];
          
          // Check if this elem was already mapped
          if ( actMap[locElem-1].GetSize() == 0 && locElem > 0 ) {
            UInt sum = 0;
            UInt max = 0;
            for( UInt iDof = 0; iDof < dofsPerElem; iDof++ ) {
              sum +=  numFcns[iDof][0];
              if( numFcns[iDof][0] > max )
                max = numFcns[iDof][0];
            }
            //LOG_DBG3(disContEqnMap) << "1-dof: " << numFcns[0].Serialize() << std::endl;
            //LOG_DBG3(disContEqnMap) << "2-dof: " << numFcns[1].Serialize() << std::endl;
            //LOG_DBG3(disContEqnMap) << "3-dof: " << numFcns[2].Serialize() << std::endl;
            //LOG_DBG3(disContEqnMap) << "maximum for this elem"  << max << std::endl;
            //LOG_DBG3(disContEqnMap) << "Elem got " << sum << " equation numbers!\n";
             

            // iterate over all element functions
            UInt pos = 0;
            for( UInt iFcn = 0; iFcn < max; iFcn++ ) {
              if( actMap[locElem-1].GetSize() == 0 ) {
                actMap[locElem-1].Resize( dofsPerElem * max );
                actMap[locElem-1].Init();
              }
              
              // iterate over all element dofs
              for( UInt iDof = 0; iDof < dofsPerElem; iDof++ ) {
                // Check if the related nodes have a 0 equations number
                // (= are (in-)hom. Dirichlet nodes)
                bool allZero = true;
                StdVector<UInt> const & elemNodes = actEl.connect;
                UInt elemIdx = 0;
                for( UInt iNode = 0; iNode < elemNodes.GetSize(); iNode++ ) {
                  for ( UInt k=0;k< mesh2DisPdeNode_[elemNodes[iNode]-1].GetSize() ;k++ ) {
                      if(static_cast<UInt>(mesh2DisPdeNode_[elemNodes[iNode]-1][k][0]) == actEl.elemNum )  
                        elemIdx=k;
                  }
                  if( actNodeMap[mesh2DisPdeNode_[elemNodes[iNode]-1][elemIdx][1]-1][iDof] != 0) {
                    allZero = false;
                    break;
                  }
                }
                if(  !allZero && iFcn < numFcns[iDof][0] ) {
                  actMap[locElem-1][pos++] = ++numEqns_;
                } else {
                  actMap[locElem-1][pos++] = 0;
                } // check for number of functions
              } // loop over dofs
            } // loop over functions
          } // check if element was already mapped
        } // loop over elements
      } // loop over entitylists
    } // loop over results
  }

  void DiscontinuousEqnMap::CalcElemConstEquations( UInt phase ) {

    // Big outer loop over all nodal mapped element lists
    ResultEntityMap::iterator listIt;
    
    // MAGIC number, which gets assignetd to all nodes,
    // which have not yet an equation number
    const Integer NO_EQN = -333;
    
    for( listIt = elemConstMappedList_.begin(); 
    listIt != elemConstMappedList_.end(); 
    listIt++ ) {

      // Remeber current result and list of elementLists
      const ResultInfo & actRes = listIt->first;
      StdVector<shared_ptr<EntityList> > & actLists = listIt->second;
      
      
      // Get grip of homogeneous and in-homogeneous boundary conditions
      // for this tpye of result
      ResultHdBcMap::iterator hdBcIt = hdBcs_.find( actRes );
      ResultIdBcMap::iterator idBcIt = idBcs_.find( actRes );
      ResultConstraintMap::iterator csIt = constraints_.find( actRes );
      
      StdVector<Vector<Integer> > & actMap = elemEqns_[actRes];

      // Get number of dofs
      UInt dofsPerElem = actRes.dofNames.GetSize();
    
      // Idea of the algorithm:
      //
      // -- PHASE 1 --
      // Step 1:  Initialize actMap with NO_EQN
      // Step 2:  For each entry in homoDirichletNodes_ set the corresponding
      //          entry in pdeNode2eqn_ to 0
      // Step 2b: For each entry in inhomoDirichletNodes_ set the corresponding
      //          entry in pdeNode2eqn_ to -1
      // Step 3:  For each entry in constraintSlaveNodes_ set the corresponding
      //          entry in pdeNode2eqn_ to 0
      // Step b:  Loop over all entries in pde2Meshnode
      //          and assign each non-zero entry an equation number

      // Step 4:  Loop over all nodes of given entity lists
      //          and assign each non-zero entry an equation number

      // Step 5:  Loop over the whole map and set all entries with NO_EQN to 0
      
      // -- PHASE 2 --
      // Step 5:  Afterwards loop again over all nodes in constraintSlaveNodes_
      //          and set the corresponding entry in pdeNode2EQN_ to the
      //          negative of the value of constraintMasterNode
      // Step 5b: Loop over all entries in inhomoDirichletNodes_ and assign that
      //          dof an equation number after the hightest equation number of
      //          the free dofs
      //
      
      if( phase == 1 ) {

        // ------
        // STEP 1
        // ------
        //UInt multipleBCs = 0;
        Matrix<UInt> countElems;
        countElems.Resize( numLocElems_, dofsPerElem );
        countElems.Init();

        actMap.Resize( numLocElems_);
        for( UInt i = 0; i < actMap.GetSize(); i++) {
          actMap[i].Resize( dofsPerElem );
          actMap[i].Init( NO_EQN );
        }

        // -------
        // STEP 2b
        // -------

        countElems.Init();

        // Check if any inhom. boundary condition is defined for the current
        // result
        if( idBcIt != idBcs_.end() ) {
          IdBcList const & actIdBcList = idBcIt->second;

          for ( UInt i = 0; i < actIdBcList.GetSize(); i++ ) {
            StdVector<UInt> nodes;
            EntityIterator elemIt = actIdBcList[i]->entities->GetIterator();

            UInt actDof = actIdBcList[i]->dof;

            for( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
              UInt actElem = elemIt.GetElem()->elemNum;
              if ( mesh2PdeElem_[ actElem - 1 ] - 1 < 0 ) {
                (*warning) << "EqnMap::CalcElemEquations: Inhom. Dirichlet "
                << "elem #" << actElem
                << " is not contained in any of the regions for "
                << "this Pde";
                Warning( __FILE__, __LINE__ );
              }
              else if ( countElems[mesh2PdeElem_[actElem-1]-1]
                                   [actDof-1] != 0 ) {
                (*warning) << "EqnMap::CalcElemEquations: Inhom. Dirichlet "
                << "elem #" << actElem
                << "\nappeared already at least once in the list of "
                << "boundary nodes for this Pde!\n Please check, if "
                << "this node is defined in more than one level of "
                << "boundary nodes!";
                Warning( __FILE__, __LINE__ );
              }
              else {
                actMap[mesh2PdeElem_[actElem-1]-1] [actDof-1] = -1;
                countElems[mesh2PdeElem_[actElem-1]-1][actDof-1]++;
                // In any case we have to increment the number of idBC-conditions
                numIdBcs_++;
              }
            }
          }
        }

        // ------
        // STEP 4
        // ------
        // Initialize countNodes to zero. It will be used to count if
        // a node got already an equation number
        countElems.Init();

        StdVector<Elem*> elems;
        // Iterate over all element list belonging to this result
        // and assign each of them a node number
        for( UInt iList = 0; iList < actLists.GetSize(); iList++ ) {

          EntityIterator it = actLists[iList]->GetIterator();

          Integer locElem = 0;
          for ( it.Begin(); !it.IsEnd(); it++ ) {
            locElem = mesh2PdeElem_[(it.GetElem()->elemNum)-1];

            // Check if this element was already mapped
            if ( actMap[locElem-1].GetSize() == 0 ) {
              actMap[locElem-1].Resize( dofsPerElem );
              actMap[locElem-1].Init();
            }

            for ( UInt iDof = 0; iDof < dofsPerElem; iDof++ ) {
              if ( actMap[locElem-1][iDof] != NO_EQN &&
                  countElems[locElem-1][iDof] == 0) {
                numEqns_++;
                actMap[locElem-1][iDof] = numEqns_;
                countElems[locElem-1][iDof] = 1;
              }
            }
          }
        }
        

      } else if( phase == 2 ) {

        // ------
        // STEP 5
        // ------
        // Check if any inhom. boundary condition is defined for the current
        // result
        Integer locElem = 0;
        if( idBcIt != idBcs_.end() ) {
          IdBcList const & actIdBcList = idBcIt->second;

          for ( UInt i = 0; i < actIdBcList.GetSize(); i++ ) {
            StdVector<UInt> nodes;
            EntityIterator elemIt = actIdBcList[i]->entities->GetIterator();

            UInt actDof = actIdBcList[i]->dof;

            for( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
              UInt actElem = elemIt.GetElem()->elemNum;
              locElem = mesh2PdeElem_[actElem-1];
              if( locElem > 0 ) { 
                if(  actMap[locElem-1] [actDof-1] == -1  ) {
                  numEqns_++;
                  actMap[locElem-1] [actDof-1] = numEqns_;
                }
                }
            }
          }
        }
      } else {
        EXCEPTION( "Phase '" << phase << "' does not exist!" );
      }
    }

  }
  
  void DiscontinuousEqnMap::CalcEdgeEquations( UInt phase ) {

    // MAGIC number, which gets assignetd to all edges,
    // which have not yet an equation number
    const Integer NO_EQN = -333;
     
    // Big outer loop over all edge mapped element lists
    ResultEntityMap::iterator listIt;

    for( listIt = edgeMappedList_.begin(); 
    listIt != edgeMappedList_.end(); 
    listIt++ ) {

      // Remeber current result and list of elementLists
      const ResultInfo & actRes = listIt->first;
      StdVector<shared_ptr<EntityList> > & actLists = listIt->second;
      StdVector<Vector<Integer> > & actMap = edgeEqns_[actRes];
      
      // Get grip of homogeneous and in-homogeneous boundary conditions
      // for this tpye of result
      ResultHdBcMap::iterator hdBcIt = hdBcs_.find( actRes );

      if( phase == 1 ) {

        UInt dofsPerEdge = actRes.dofNames.GetSize();
        actMap.Resize( numLocEdges_  );

        // Fetch also nodal map
        Matrix<Integer> & actNodeMap = nodeEqns_[actRes];

        // --------
        //  Step 1
        // --------
        if( hdBcIt != hdBcs_.end() ) {
          HdBcList const & actHdBcList = hdBcIt->second;
          Integer locEdge = NO_EQN;
          for ( UInt i = 0; i < actHdBcList.GetSize(); i++ ) {

            // check, if entitylist consists of elements
            if( actHdBcList[i]->entities->GetType() == EntityList::ELEM_LIST ||
                actHdBcList[i]->entities->GetType() == EntityList::SURF_ELEM_LIST ) {
              EntityIterator it = actHdBcList[i]->entities->GetIterator();

              // iterate over all elements of this list
              for( it.Begin(); !it.IsEnd(); it++ ) {

                // obtain edges
                StdVector<Integer> const & edges = it.GetElem()->edges;
                for( UInt iEdge = 0; iEdge < edges.GetSize(); iEdge++ ) {

                  locEdge = NO_EQN;
                  for (UInt k = 0 ;k < mesh2DisPdeEdge_[ std::abs(edges[iEdge])-1 ].GetSize() ; k++)
                  {
                    if(static_cast<UInt>(mesh2DisPdeEdge_[ std::abs(edges[iEdge])-1 ][k][0]) == it.GetElem()->elemNum)
                      locEdge = mesh2DisPdeEdge_[std::abs(edges[iEdge]) - 1][k][1];
                  }

                  if( locEdge > 0 ) {
                    actMap[locEdge-1].Resize(dofsPerEdge);
                    actMap[locEdge-1].Init();
                  }
                }
              }
            }
          }
        }

        // Iterate over all entitylists
        for( UInt iList = 0; iList < actLists.GetSize(); iList++ ) {
          EntityIterator it = actLists[iList]->GetIterator();
  
          // Iterate over all elements within this list
          for( it.Begin(); !it.IsEnd(); it++ ) {

            // Get grip of element
            const Elem & actEl = *(it.GetElem());
                        
            // Get number of unknowns for each edge
            StdVector< StdVector<UInt> > numFcns; 
            numFcns.Resize( dofsPerEdge );

            // iterate over all dofs of this result
            for( UInt iDof = 0; iDof < dofsPerEdge; iDof++ ) {
              actEl.ptElem->GetNumFncs( numFcns[iDof], actRes.fctType, 
                                        AnsatzFct::EDGE, iDof );
              //LOG_DBG3(disContEqnMap) << "numFncs, dof " << iDof+1 << ": " 
              //          << numFcns[iDof].Serialize() << std::endl;
              // assert that we have as many entries as we have edges
              assert( numFcns[iDof].GetSize() == actEl.edges.GetSize() );
            }
            Integer locEdge = -999;
            // Iterate over all edges of this element
            for( UInt iEdge = 0; iEdge < actEl.edges.GetSize(); iEdge++ ) 
            {
              // determine local edge
              locEdge = -999;
              for ( UInt k = 0;k< mesh2DisPdeEdge_[std::abs(actEl.edges[iEdge])-1].GetSize() ;k++ ) {
                if(static_cast<UInt>(mesh2DisPdeEdge_[std::abs(actEl.edges[iEdge]) - 1][k][0]) == actEl.elemNum ) {
                  locEdge = mesh2DisPdeEdge_[std::abs(actEl.edges[iEdge]) - 1][k][1];
                }
              }
              if(locEdge == -999){
                LOG_DBG3(disContEqnMap) << "local Edge not found, check the mesh2DisPdeEdge_ Array" << std::endl;
              }

              Edge const & edge = ptGrid_->GetEdge( std::abs(actEl.edges[iEdge] ));
              
              // Check if this edge was already mapped
              if ( actMap[locEdge-1].GetSize() == 0 && locEdge > 0 ) {
                
                // sum up unknowns of this dof
                UInt sum = 0;
                UInt max = 0;
                for( UInt iDof = 0; iDof < dofsPerEdge; iDof++ ) {
                  sum +=  numFcns[iDof][iEdge];
                  if( numFcns[iDof][iEdge] > max )
                    max = numFcns[iDof][iEdge];
                }
                //LOG_DBG3(disContEqnMap) << "1-dof: " << numFcns[0].Serialize();
                //LOG_DBG3(disContEqnMap) << "2-dof: " << numFcns[1].Serialize();
                //LOG_DBG3(disContEqnMap) << "3-dof: " << numFcns[2].Serialize();
                LOG_DBG3(disContEqnMap) << "maximum for edge " << iEdge 
                                 << ": " << max;
                LOG_DBG3(disContEqnMap) << "Edge got " << sum 
                                 << " equation numbers!" << std::endl;
                
                // iterate over all edge functions
                UInt pos = 0;
                UInt counter = 0;
                for( UInt iFcn = 0; iFcn < max; iFcn++ ) {
                  if( actMap[locEdge-1].GetSize() == 0 ) {
                    actMap[locEdge-1].Resize( dofsPerEdge * max );
                    actMap[locEdge-1].Init();
                  }
                  
                  Integer elemIdx1 = 0;
                  Integer elemIdx2 = 0;
                  bool foundIdx1 = false;
                  bool foundIdx2 = false;
                  // iterate over all dofs
                  for( UInt iDof = 0; iDof < dofsPerEdge; iDof++ ) {
                    elemIdx1 = -333;
                    elemIdx2 = -333;
                    foundIdx1 = false;
                    foundIdx2 = false;
                    //find the correct element the edge belons to
                    for ( UInt k=0;k< mesh2DisPdeNode_[edge.nodes[0]-1].GetSize() ;k++ ) {
                      if(static_cast<UInt>(mesh2DisPdeNode_[edge.nodes[0]-1][k][0])== actEl.elemNum ) {
                        elemIdx1=k;
                        foundIdx1 = true;
                      }
                    }
                    for ( UInt k=0;k< mesh2DisPdeNode_[edge.nodes[1]-1].GetSize() ;k++ ) {
                      if(static_cast<UInt>(mesh2DisPdeNode_[edge.nodes[1]-1][k][0])== actEl.elemNum ) {
                        elemIdx2=k;
                        foundIdx2 = true;
                      }
                    }
                    assert(foundIdx1); 
                    assert(foundIdx2); 

                    //spectral and legendre functions have to be treated
                    //differently in the case of inhomognious BCs due to the
                    //non-hirarchical structure of sprectral approximations
                    if( actRes.fctType->GetType() == AnsatzFct::SPECTRAL) {
                      // Check if the related nodes have a 0 equations number
                      // (= are (in-)hom. Dirichlet nodes)
                      // of if the edge with this dof has a smaller order
                      // than the maximum for this edge -> assign 0 equation number
                      if(  (actNodeMap[mesh2DisPdeNode_[edge.nodes[0]-1][elemIdx1][1]-1][iDof] != 0  || 
                            actNodeMap[mesh2DisPdeNode_[edge.nodes[1]-1][elemIdx2][1]-1][iDof] != 0) &&
                           (iFcn < numFcns[iDof][iEdge]) ) 
                      {
                        actMap[locEdge-1][pos++] = numEqns_ + counter + 1;
                        if( actNodeMap[mesh2DisPdeNode_[edge.nodes[0]-1][elemIdx1][1]-1][iDof] < 0  && 
                            actNodeMap[mesh2DisPdeNode_[edge.nodes[1]-1][elemIdx2][1]-1][iDof] < 0)
                          numIdBcs_++;

                        counter++;
                      }else
                        actMap[locEdge-1][pos++] = 0;
                    } else if ( actRes.fctType->GetType() == AnsatzFct::LEGENDRE ) {
                      if(  (actNodeMap[mesh2DisPdeNode_[edge.nodes[0]-1][elemIdx1][1]-1][iDof] > 0 || 
                            actNodeMap[mesh2DisPdeNode_[edge.nodes[0]-1][elemIdx2][1]-1][iDof] > 0 ) && 
                           (iFcn < numFcns[iDof][iEdge]) ) {
                        actMap[locEdge-1][pos++] = ++numEqns_;
                      } else {
                        actMap[locEdge-1][pos++] = 0;
                      }// check if edge has one non-zero node
                    } else {
                      // This here is the Nedelec case -> just number the edges
                      actMap[locEdge-1][pos++] = ++numEqns_;
                    }
                  } // loop over dofs
                } // loop over functions
                numEqns_ += counter;
                counter = 0;
              } // check if edge was already mapped
            } // loop over all edges
          } // looop over all elements
        } // loop over all entitylists

      } else {
        // === PHASE 2 ===  
      } 
    } // loop over all results
  } 

  void DiscontinuousEqnMap::CalcFaceEquations( UInt phase ) {

    LOG_TRACE(disContEqnMap) << "Starting to map face equations\n";

    // Big outer loop over all face mapped element lists
    ResultEntityMap::iterator listIt;

    for( listIt = faceMappedList_.begin(); 
         listIt != faceMappedList_.end(); 
         listIt++ ) {

      // Remeber current result and list of elementLists
      const ResultInfo & actRes = listIt->first;
      StdVector<shared_ptr<EntityList> > & actLists = listIt->second;
      StdVector<Vector<Integer> > & actMap = faceEqns_[actRes];
      
      UInt dofsPerFace = actRes.dofNames.GetSize();
      actMap.Resize( numLocFaces_  );

      // Fetch also nodal map
      Matrix<Integer> & actNodeMap = nodeEqns_[actRes];

      Integer locFace = 0;

      // Iterate over all entitylists
      for( UInt iList = 0; iList < actLists.GetSize(); iList++ ) {

        EntityIterator it = actLists[iList]->GetIterator();

        // Iterate over all elements within this list
        for( it.Begin(); !it.IsEnd(); it++ ) {

          // Get grip of element
          const Elem & actEl = *(it.GetElem());

          // Get number of unknowns for each face
          StdVector< StdVector<UInt> > numFcns; 
          numFcns.Resize( dofsPerFace );

          // iterate over all dofs of this result
          for( UInt iDof = 0; iDof < dofsPerFace; iDof++ ) {
            actEl.ptElem->GetNumFncs( numFcns[iDof], actRes.fctType, 
                                      AnsatzFct::FACE, iDof );
          }
          
          // Iterate over all faces of this element
          for( UInt iFace = 0; iFace < actEl.faces.GetSize(); iFace++ ) {
            
            // Check if there are any faces at all present
            //LOG_DBG3(disContEqnMap) << "Element has " << numFcns[iFace]
            //                 << "face functions!";
            //if( numFcns[iFace].GetSize() == 0 ) {
            //  continue;
            //}
            
            // determine local face number
            UInt locElemPos = 999;
            locFace = -999;
            for ( UInt k = 0;k< mesh2DisPdeFace_[std::abs(actEl.faces[iFace])-1].GetSize() ;k++ )
            {
              if(static_cast<UInt>(mesh2DisPdeFace_[std::abs(actEl.faces[iFace]) - 1][k][0]) == actEl.elemNum )
              {
                locFace = mesh2DisPdeFace_[std::abs(actEl.faces[iFace]) - 1][k][1];
                locElemPos = k;
              }
            }
            if(locFace == -999 || locElemPos == 999){
              LOG_DBG3(disContEqnMap) << "local Edge not found, check the mesh2DisPdeEdge_ Array" << std::endl;
            }

            Face const & face = ptGrid_->GetFace( actEl.faces[iFace] );
            LOG_DBG3(disContEqnMap) << "Actual face: " << actEl.faces[iFace] << " / "
                             << locFace << " (global / local )";
            
            // Check if this face was already mapped
            if ( actMap[locFace-1].GetSize() == 0 && locFace > 0 ) {
              
              // sum up unknowns of this dof
              UInt sum = 0;
              UInt max = 0;
              for( UInt iDof = 0; iDof < dofsPerFace; iDof++ ) {
                sum +=  numFcns[iDof][iFace];
                if( numFcns[iDof][iFace] > max )
                  max = numFcns[iDof][iFace];
              }
              //LOG_DBG3(disContEqnMap) << "1-dof: " << numFcns[0].Serialize();
              //LOG_DBG3(disContEqnMap) << "2-dof: " << numFcns[1].Serialize();
              //LOG_DBG3(disContEqnMap) << "3-dof: " << numFcns[2].Serialize();
              LOG_DBG3(disContEqnMap) << "maximum for face " << iFace 
                               << ": " << max;
              LOG_DBG3(disContEqnMap) << "Face got " << sum << " equation numbers!\n";
             
              
              // iterate over all functions of this face
              UInt pos = 0;
              for( UInt iFcn = 0; iFcn < max; iFcn++ ) {
                
                if( actMap[locFace-1].GetSize() == 0 ) {
                  actMap[locFace-1].Resize( dofsPerFace * max);
                  actMap[locFace-1].Init();
                }
              
                // iterate over all dofs
                for( UInt iDof = 0; iDof < dofsPerFace; iDof++ ) {
                  
                  // Check if the related nodes have an equation number =< 0
                  // (= are (in-)hom. Dirichlet nodes)
                  bool allFixed = true;
                  bool isInHomBoundaryFace = true;
                  if( actRes.fctType->GetType() == AnsatzFct::SPECTRAL) {
                    for( UInt iNode = 0; iNode < face.nodes.GetSize(); iNode++ ) {
                      for ( UInt k=0;k< mesh2DisPdeNode_[face.nodes[iNode]-1].GetSize() ;k++ ) {
                        if(static_cast<UInt>(mesh2DisPdeNode_[face.nodes[iNode]-1][k][0])== actEl.elemNum ){
                          locElemPos=k;
                          break;
                        }
                      }
                      if( actNodeMap[mesh2DisPdeNode_[face.nodes[iNode]-1][locElemPos][1]-1][iDof] != 0) {
                        allFixed = false;
                      }
                      if( actNodeMap[mesh2DisPdeNode_[face.nodes[iNode]-1][locElemPos][1]-1][iDof] > 0) {
                        isInHomBoundaryFace = false;
                      }
                    }
                    if( !allFixed  && (iFcn < numFcns[iDof][iFace]) ) 
                    {
                      actMap[locFace-1][pos++] = ++numEqns_;
                      LOG_DBG3(disContEqnMap) << "Face #" << actEl.faces[iFace]
                                       << " got equation number " << numEqns_-1 << "For Discontinuous Element #" << actEl.elemNum;
                      if( isInHomBoundaryFace )
                        numIdBcs_++;

                    }else{
                      actMap[locFace-1][pos++] = 0;
                      LOG_DBG3(disContEqnMap) << "-> Nodal equations 0 for dof" << iDof+1;
                    } // check for zero nodes
                  } else{
                    for( UInt iNode = 0; iNode < face.nodes.GetSize(); iNode++ ) {
                      for ( UInt k=0;k< mesh2DisPdeNode_[face.nodes[iNode]-1].GetSize() ;k++ ) {
                        if(static_cast<UInt>(mesh2DisPdeNode_[face.nodes[iNode]-1][k][0])== actEl.elemNum ){
                          locElemPos=k;
                          break;
                        }
                      }
                      if( actNodeMap[mesh2DisPdeNode_[face.nodes[iNode]-1][locElemPos][1]-1][iDof] > 0) {
                        allFixed = false;
                        break;
                      }
                    }
                    if( !allFixed && (iFcn < numFcns[iDof][iFace]) ) {
                      actMap[locFace-1][pos++] = ++numEqns_;
                      LOG_DBG3(disContEqnMap) << "Face #" << actEl.faces[iFace]
                                       << " got equation number " << numEqns_-1;
                    } else {
                      actMap[locFace-1][pos++] = 0;
                      LOG_DBG3(disContEqnMap) << "-> Nodal equations 0 for dof" << iDof+1;
                    } // check for zero nodes
                  }// spectral or not
                } // loop over dofs
              }  // loop over functions
            } else {
              LOG_DBG3(disContEqnMap) << "--> Already mapped\n ";
            }// check if already mapped
          } // loop over faces
        } // loop over elements
      } // loop over results
    }
    
    LOG_TRACE(disContEqnMap) << "Finished mapping face equations\n";
    
  }

  void DiscontinuousEqnMap::GetNodesOfEntities( StdVector<UInt>& nodes,
      shared_ptr<EntityList> ent ) {
    
    // Get type of entries of the particular entity list
    //EntityList::ListType type = ent->GetType();
    
    shared_ptr<ElemList> elemList;
    shared_ptr<SurfElemList> sElemList;
    shared_ptr<NodeList> nodeList;
    shared_ptr<RegionList> regionList;
    StdVector<UInt> helpNodes;
    EntityIterator it;
    
    // get name of entitylist
    std::string name= ent->GetName();
    ptGrid_->GetNodesByName( nodes, name );

    //    switch( type ) {
    //      
    //    case EntityList::ELEM_LIST:
    //      elemList= 
    //  dynamic_pointer_cast<ElemList, EntityList>(ent);
    //      ptGrid_->GetNodesByRegion( nodes, elemList->GetRegion() );
    //      break;
    //      
    //    case EntityList::SURF_ELEM_LIST:
    //      sElemList = 
    //  dynamic_pointer_cast<SurfElemList, EntityList>(ent);
    //      ptGrid_->GetNodesByRegion( nodes, sElemList->GetRegion() );
    //      break;
    //
    //    case EntityList::REGION_LIST:
    //      regionList= 
    //  dynamic_pointer_cast<RegionList, EntityList>(ent);
    //      it = regionList->GetIterator();
    //      for( ; !it.IsEnd(); it++ ) {
    //        helpNodes.Clear();
    //        ptGrid_->GetNodesByRegion( helpNodes, it.GetRegion() );
    //        for( UInt i = 0; i < helpNodes.GetSize(); i++ ) {
    //          nodes.Push_back( helpNodes[i] );
    //        }
    //      }
    //      break;
    //
    //    case EntityList::NODE_LIST:
    //      nodeList = 
    //  dynamic_pointer_cast<NodeList, EntityList>(ent);        
    //      nodes = nodeList->GetNodes();
    //      break;
    //      
    //    default :
    //      std::string listString;
    //      EntityList::Enum2String( type, listString );
    //      EXCEPTION( "'" << listString
    //                 << "' is no EntityList with nodal information." );
    //      
    //    }
  }

  /*void DiscontinuousEqnMap::GetResEqns(  StdVector<Integer>& eqns, const ResultInfo& result ) const 
  {
    StdVector< shared_ptr< EntityList > > eLists = (resEntMap_.find(result) )->second;
    StdVector< Integer > actEqns;
    for ( UInt i= 0;i < eLists.GetSize() ;i++ )
    {
      EntityIterator entIt = eLists[i]->GetIterator();
      UInt size = eLists[i]->GetSize();
      entIt.Begin();
      while ( !entIt.IsEnd() )
      {
        GetEqns(actEqns,result,entIt);
        for ( UInt x = 0; x<actEqns.GetSize() ; x++ )
        {
          if(eqns.Find(actEqns[x]) == -1)
            eqns.Push_back(actEqns[x]);
        }
        entIt++;
      }
    }
  }*/
}
