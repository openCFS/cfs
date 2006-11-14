#include "nodestoresol.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField {


  // ***********************
  //   Default Constructor
  // ***********************
  template<class TYPE>
  NodeStoreSol<TYPE>::NodeStoreSol() {

    ENTER_FCN( "NodeStoreSol::NodeStoreSol()" );
  
    numNodes_     = 0;
    numSolutions_ = 0;
    length_       = 0;
    totalDofs_    = 0;
  }


  // ***************************
  //   Constructor (Version 1)
  // ***************************
  template<class TYPE>
  NodeStoreSol<TYPE>::NodeStoreSol( UInt numNodes, 
                                    StdVector<SolutionType> solTypes, 
                                    StdVector<UInt> solDofs ) {
    ENTER_FCN( "NodeStoreSol::NodeStoreSol(numNodes, solTypes, solDofs" );
    Error( "Not implemented here", __FILE__, __LINE__ );
  }


  // ***************************
  //   Constructor (Version 2)
  // ***************************
  template<class TYPE>
  NodeStoreSol<TYPE>::NodeStoreSol( const UInt numNodes,
                                    const SolutionType solType,
                                    const UInt numDofs ) {
    ENTER_FCN( "NodeStoreSol::NodeStoreSol(numNodes, solType, soDofs" );
    Error( "Not implemented here", __FILE__, __LINE__ );
  }


  // ***************************
  //   Constructor (Version 3)
  // ***************************
  template<class TYPE>
  NodeStoreSol<TYPE>::NodeStoreSol( const NodeStoreSol &x ) {
    ENTER_FCN( "NodeStoreSol::NodeStoreSol(const NodeStoreSol)" );
    Error( "Not implemented here", __FILE__, __LINE__ );
  }


  // **************
  //   Destructor
  // **************
  template<class TYPE>
  NodeStoreSol<TYPE>::~NodeStoreSol() {
    ENTER_FCN( "NodeStoreSol::~NodeStoreSol" );
  }


  // *************
  //   IsComplex
  // *************
  template<class TYPE>
  bool NodeStoreSol<TYPE>::IsComplex() {
    return false;
  }

  template<>
  bool NodeStoreSol<Complex>::IsComplex() {
    return true;
  }


  // *********
  //   Clear
  // *********
  template<class TYPE>
  void NodeStoreSol<TYPE>::Clear() {
    ENTER_FCN( "NodeStoreSol::Clear" );
    Error( "NodeStoreSol::Clear() not implemented", __FILE__, __LINE__ );
  }


  // ********
  //   Init
  // ********
  template<class TYPE>
  void NodeStoreSol<TYPE>::Init() {
    ENTER_FCN( "NodeStoreSol::Init" );
    Init( TYPE() );
  }


  // *****************
  //   SetPtrEQNData
  // *****************
  template<class TYPE>
  void NodeStoreSol<TYPE>::SetPtrEQNData( EqnMap * eqnMap,
                                          Grid *ptGrid ) {
    ENTER_FCN( "NodeStoreSol::SetPtrEQNData" );
    ptGrid_ = ptGrid;
    eqnMap_ = eqnMap;
  }


  // ********
  //   Init
  // ********
  template<class TYPE>
  void NodeStoreSol<TYPE>::Init( const TYPE val ) {

    ENTER_FCN( "NodeStoreSol::Init" );
  
#ifdef CHECK_INITIALIZED
    if (numSolutions_ == 0 || numNodes_ == 0 || solTypes_.size() == 0 ||
        solDofs_.size() == 0 ) {
      (*error) << "NodeStoreSol::Init:"
               << "\n numSolutions ....... " << numSolutions_
               << "\n numNodes ........... " << numNodes_
               << "\n totalDofs .......... " << totalDofs_
               << "\n length ............. " << length_
               << "\n size of solDofs .... " << solDofs_.size()
               << "\n size of solTypes ... " << solTypes_.size()
               << "\n size of offsets .... " << solOffset_.size()
               << "\n Before calling Init(), the number of solutions, "
               << "nodes, types and dofs has to be set to NON-ZERO value!";
      Error( __FILE__, __LINE__ );
    }
#endif
  
    // only the first time the struct gets initialized
    if ( length_ == 0 ) {
      if ( solDofs_.size() != numSolutions_ ||
           solTypes_.size() != numSolutions_ ) {
        (*error) << "Inconsistent definition of Storesolution class. "
                 << "Maybe you have to call 'Clear()' before using a "
                 << "modified data layout!";
        Error( __FILE__, __LINE__ );
      }

      totalDofs_ = 0;
      
      // iterate over all solutiontypes and
      // sum up their number of dof
      std::map<SolutionType,UInt>::iterator it;
      for ( it = solDofs_.begin(); it != solDofs_.end(); it++ ) {
        // set offset of current solution
        // w.r.t. to starting position
        solOffset_[(*it).first] = totalDofs_;
        totalDofs_ += (*it).second;     
      }

      length_ = totalDofs_ * numNodes_;

      // Check for NULL-Pointer
      if ( eqnMap_ == NULL ) {
        Error("NodeStoreSol::Init: Pointer to EQN is NULL!",
              __FILE__, __LINE__);
      }
      
      // Determine the length of the solution vector
      lengthVector_ = eqnMap_->GetNumEqns();
      eqnDofs_ = 1;

      data_.Resize(lengthVector_);
      data_.Init(val);
    }
  }


  // *******************
  //   SetNumSolutions
  // *******************
  template<class TYPE>
  void NodeStoreSol<TYPE>::SetNumSolutions( const UInt nSols ) {
    ENTER_IFCN( "NodeStoreSol::SetNumSolutions" );
    numSolutions_ = nSols;
    length_ = 0;
  }


  // ***************
  //   SetNumNodes
  // ***************
  template<class TYPE>
  void NodeStoreSol<TYPE>::SetNumNodes( const UInt nNodes ) {
    ENTER_IFCN( "NodeStoreSol::SetNumNodes" );
    numNodes_ = nNodes;
    length_ = 0;
  }


  // *******************
  //   SetSolutionType
  // *******************
  template<class TYPE>
  void NodeStoreSol<TYPE>::SetSolutionType( const SolutionType solType,
                                            const UInt numSolution ) {

    ENTER_IFCN( "NodeStoreSol::SetSolutionType" );

#ifdef CHECK_INDEX
    if (numSolution >= numSolutions_) {
      Error("NodeStoreSol: Index out of Bounds",__FILE__,__LINE__);
    }
#endif

    // Check, if the object contains only one entry and
    // if this only entry is overwritten
    if ( solTypes_.size() == 1 && numSolutions_ == 1 && numSolution == 0 ) {
      solTypes_.clear();
    }

    solTypes_[solType] = numSolution;
    length_ = 0;
  }


  // **************
  //   SetNumDofs
  // **************
  template<class TYPE>
  void NodeStoreSol<TYPE>::SetNumDofs( const UInt dof,
                                       const SolutionType sol ) {

    ENTER_IFCN( "NodeStoreSol::SetDof" );

    // check if only one dof was assigned
    // -> only one entry exists
    if ( numSolutions_ == 1 && sol == NO_SOLUTION_TYPE ) {
      solDofs_.clear();
      solDofs_[(*solTypes_.begin()).first] = dof;
    }
    else {
      solDofs_[sol] = dof;
    }
    length_ = 0;
  }


  // *******************
  //   GetSolutionType
  // *******************
  template<class TYPE>
  void NodeStoreSol<TYPE>::
  GetSolutionTypes( StdVector<SolutionType> &solTypes ) const {

    ENTER_FCN( "NodeStoreSol::GetSolutionTypes" );

    solTypes.Resize(numSolutions_);
    std::map<SolutionType,UInt>::const_iterator it;
 
    for ( it = solTypes_.begin(); it != solTypes_.end(); it++ ) {
      solTypes[(*it).second] = (*it).first;
    }
  }

  template<class TYPE>
  TYPE& NodeStoreSol<TYPE>::operator()(UInt node, UInt dof)
  {
    ENTER_IFCN("NodeStoreSol::operator()");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
    if (numSolutions_ > 1) 
      Error("NodeStoreSol:operator(): Only defined objects with one type of solution!",__FILE__,__LINE__);
#endif
    Integer eqnNr;
    eqnNr = eqnMap_->GetNodeEqn( node+1, dof+1 );

    if (eqnNr > 0)
      return data_[abs(eqnNr)-1];
    else
      Error("NodeStoreSol::operator(): This operator gives only writing access to non-BC nodes.",
            __FILE__, __LINE__);
  }

  template<class TYPE>
  TYPE  NodeStoreSol<TYPE>:: operator()(UInt node, UInt dof) const
  {
    ENTER_IFCN("NodeStoreSol::operator()");
  
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
    if (numSolutions_ > 1) 
      Error("NodeStoreSol:operator(): Only defined objects with one type of solution!",__FILE__,__LINE__);
#endif
  
    Integer eqnNr;
    eqnNr = eqnMap_->GetNodeEqn( node+1,dof+1 );

    Warning ("Is this operator ever be used?", __FILE__, __LINE__);
    if (eqnNr > 0)
      return data_[abs(eqnNr)-1];
    else
      return TYPE();
  }


  template<class TYPE>
  void NodeStoreSol<TYPE>::GetGlobalSolVector(const SolutionType type,
                                              CFSVector & val) const
  {
    ENTER_FCN("NodeStoreSol::GetGlobalSolVector");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) 
      Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif


    UInt dof = (*solDofs_.find(type)).second;

    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);

    UInt globNumNodes = ptGrid_->GetNumNodes();

    temp.Resize(globNumNodes *dof);
    temp.Init();

    UInt numDofs = result_->dofNames.GetSize();
    StdVector<Integer> eqns;
    
    // Iterate over all entities
    for( UInt i = 0; i < elems_.GetSize(); i++ ) {
      
      // Get iterator
      EntityIterator it = elems_[i]->GetIterator();

      for (it.Begin(); !it.IsEnd(); it++ ) {

        // Get connectivity
        StdVector<UInt> const & nodes  = it.GetElem()->connect;

        // Get related equation numbers
        //eqnMap_->GetEqns( eqns, *result_, it );
        eqnMap_->GetNodeEqn( nodes, eqns );

        // Iterate over all nodes
        for( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
          Integer globNode = nodes[iNode];
          // Iterate over all dofs
          for (UInt iDof = 0; iDof < numDofs; iDof++ ) {
            
            Integer eqnNr = eqns[iNode*numDofs+iDof];

            if( eqnNr != 0 ) {
              temp[(globNode-1)*numDofs + iDof] = data_[std::abs(eqnNr-1)];
            } else if ( eqnNr == 0 ) {
              temp[(globNode-1)*numDofs + iDof] = TYPE();
            }
          }
        }
        

      }
    }
    

    // OLD-version
//     dofsPerEQN = ptEQN_->GetNumDofsPerEQN();

//     UInt numLocNodes = ptEQN_->GetNumLocalNodes();
//     UInt globNode = 0;
//     for (UInt iNode=0; iNode<numLocNodes; iNode++) {
//       for (UInt iDof=0; iDof<dof; iDof++) {
//         globNode = ptEQN_->PDE2MeshNode(iNode+1);
//         ptEQN_->Node2EQN(globNode,iDof+1+offset,eqnNr,eqnDof);
        
//         if (eqnNr > 0) 
//           temp[(globNode-1)*dof + iDof] = data_[(eqnNr-1) * dofsPerEQN + (eqnDof-1)];
//         else if (eqnNr == 0)
//           temp[(globNode-1)*dof + iDof] = TYPE();
//         else
//           temp[(globNode-1)*dof + iDof] = data_[(abs(eqnNr)-1) * dofsPerEQN + (eqnDof-1)];
        
//       }
//     } 
  }


  template<class TYPE>
  void NodeStoreSol<TYPE>::SetNodalResult(const UInt nodeNr, const CFSVector &val)
  {
    ENTER_FCN("NodeStoreSol::SetNodalResult");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
#ifdef CHECK_INDEX
    if (nodeNr >= numNodes_) 
      Error("NodeStoreSol::SetNodalResult(): index out of bounds",__FILE__,__LINE__);
    if (val.GetSize() != totalDofs_)
      Error("NodeStoreSol::SetNodalResult(): vector of incompatible dimension",__FILE__,__LINE__);
#endif

    const Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
    for (UInt i=0; i<temp.GetSize(); i++) {
      UInt eqnNr = eqnMap_->GetNodeEqn( nodeNr, i+1 );
      data_[eqnNr-1] = temp[i];
    }
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::GetNodalResult(const UInt nodeNr, CFSVector &val) const
  {
    ENTER_FCN("NodeStoreSol::GetNodalResult");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
    Error("Not implemented here", __FILE__,__LINE__);
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::GetGlobalSolVectorSingleDof(const SolutionType type, 
                                                       const UInt dof, 
                                                       CFSVector & val) const
  {
    ENTER_FCN("NodeStoreSol::GetSolVectorSingleDof");
    UInt solDof = (*solDofs_.find(type)).second;
#ifdef CHECK_INITIALIZED

    if (length_ == 0) 
      Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);


    std::string errMsg;
    if (dof > solDof)
      {
        errMsg = "NodeStoreSol::GetSolVectorSingleDof: Desired dof ";
        errMsg += GenStr(dof);
        errMsg += " is higher than dofs of solution = ";
        errMsg += GenStr(solDof);
        Error (errMsg.c_str(), __FILE__, __LINE__);
      }
#endif
  
    UInt offset = (*solOffset_.find(type)).second;
    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
    
    Integer eqnNr;
    UInt globNumNodes = ptGrid_->GetNumNodes();
    temp.Resize(globNumNodes);
    temp.Init();
    UInt numLocNodes = eqnMap_->GetNumLocalNodes();
    UInt globNode = 0;
    for (UInt iNode=0; iNode<numLocNodes; iNode++) {
      globNode = eqnMap_->Pde2MeshNode(iNode+1);
      eqnNr = eqnMap_->GetNodeEqn( globNode, dof+1+offset ); 
 
      if (eqnNr != 0) 
        temp[globNode-1] = data_[(eqnNr-1)];
      else if (eqnNr == 0)
        temp[globNode-1] = TYPE();
      else
        Error("Constraints not yet implemented!",__FILE__,__LINE__);
    } 
  
  }  


  template<class TYPE>
  void NodeStoreSol<TYPE>::GetGlobalSolVectorSingleDof(const UInt dof, 
                                                       CFSVector & val) const{
    ENTER_FCN("NodeStoreSol::GetSolVectorSingleDof");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) {
      Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
    }
#endif
  
    //Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
 
    Error ("Not working properly at the moment", __FILE__, __LINE__);

    //   UInt globalPos, eqnDofs;
  
    //   temp.Resize(ptEQN_->GetNumGlobalNodes());
    //   eqnDofs = ptEQN_->GetNumDofs();
 
    //   // Loop over all Equations
    //   for (UInt iEQN=0; iEQN<ptEQN_->GetNumEQNs(); iEQN++)
    //     {
    //       ptEQN_->EQN2SolVectorPos(iEQN+1, 1, globalPos);

    //       // ONLY TEMPORARY, until superblocknumbering
    //       // for piezoPDE works ;-)
    //       globalPos =(UInt) globalPos/eqnDofs;
      
    //       temp[globalPos] = data_[iEQN*totalDofs_ + dof];
    //     }
  }



  // *******
  //   Get
  // *******
  template<class TYPE>
  void NodeStoreSol<TYPE>::Get( const UInt nodeNr, const UInt dof,
                                TYPE &ret ) const {

    ENTER_FCN( "NodeStoreSol::Get" );

#ifdef CHECK_INITIALIZED
    if (length_ == 0) {
      Error( "NodeStoreSol::Get(): Use of uninitialized object!", __FILE__,
             __LINE__ );
    }
#endif

#ifdef CHECK_INDEX
    if ( numSolutions_ > 1 ) {
      Error( "NodeStoreSol::Get(): Only used for single solution objects!",
             __FILE__, __LINE__ );
    }
#endif

    Integer eqnNr;
    eqnNr = eqnMap_->GetNodeEqn( nodeNr + 1, dof + 1 );
    if ( eqnNr != 0 ) {
      ret = data_[ abs(eqnNr)- 1];
    }
    else {
      ret =  TYPE();
    }
  }


  // *******
  //   Get
  // *******
  template<class TYPE>
  void NodeStoreSol<TYPE>::Get(const SolutionType type, 
                               const UInt nodeNr, 
                               const UInt dof, 
                               TYPE & ret) const
  {
    ENTER_FCN("NodeStoreSol::Get");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

    Integer eqnNr;
  
    eqnNr = eqnMap_->GetNodeEqn( *result_, nodeNr+1, dof+1);

    if (eqnNr != 0)
      ret = data_[abs(eqnNr)-1];
    else
      ret =  TYPE();
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::Set(const SolutionType type, const UInt nodeNr, const UInt dof, const TYPE val)
  {
    ENTER_FCN("NodeStoreSol::Set");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
    //    Error("Not implemented here", __FILE__,__LINE__);

    UInt offset = (*solOffset_.find(type)).second;

    Integer eqnNr;
  
    eqnNr = eqnMap_->GetNodeEqn( nodeNr+1, offset+dof+1 );

    if (eqnNr != 0)
      data_[abs(eqnNr-1)] = val;
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::Add(const SolutionType type, const UInt nodeNr, const UInt dof, const TYPE val) const
  {
    ENTER_FCN("NodeStoreSol::Set");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
    Error("Not implemented here", __FILE__,__LINE__);
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::SetAlgSysVector(const CFSVector & val)
  {
    ENTER_FCN("NodeStoreSol::SetAlgSysVector");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
#ifdef CHECK_INDEX
   //  if (val.GetSize() !=  lengthVector_)
//       {
//         std::cerr << "Vector has Size" << lengthVector_ << std::endl;
//         std::cerr << "Your vector has size " << val.GetSize() << std::endl;
//         Error("NodeStoreSol::SetAlgSysVector(): Vector has wrong size!",__FILE__,__LINE__);
//       }
#endif

    const  Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
    data_ = temp;
  }
  
  template<class TYPE> 
  void NodeStoreSol<TYPE>::GetAlgSysVector(CFSVector & val) const
  {
    ENTER_FCN("NodeStoreSol::GetAlgSysVector");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
    temp = data_;
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::GetAlgSysVectorPointer(CFSVector* &ptrToVec)
  {
    ENTER_FCN("NodeStoreSol::GetAlgSysVectorPointer");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

    ptrToVec = (CFSVector*) &data_;
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::GetAlgSysVectorPointer(Vector<TYPE>* &ptrToVec)
  {
    ENTER_FCN("NodeStoreSol::GetAlgSysVectorPointer");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

    ptrToVec =  &data_;
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::CopyFromAlgSysDataPointer(Double * ptr)
  {
    ENTER_FCN("NodeStoreSol::CopyFromAlgSysDataPointer");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
    for (UInt i=0; i<lengthVector_; i++)
      {
        data_[i] = (TYPE) ptr[i];
      }

  }

  // Specialisation for complex values
  template<>
  void NodeStoreSol<Complex>::CopyFromAlgSysDataPointer(Double * ptr)
  {
    ENTER_FCN("NodeStoreSol::CopyFromAlgSysDataPointer");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
    for (UInt i=0; i<lengthVector_; i++)
      {
        data_[i] = Complex(ptr[2*i],ptr[2*i+1]);
      }

  }


  template<class TYPE>
  void NodeStoreSol<TYPE>::SetAlgSysDataPointer( UInt size, TYPE * ptr)
  {
    ENTER_FCN("NodeStoreSol::SetAlgSysDataPointer");
  
    data_.Replace( size, ptr, false );

  }


  template<>
  Double* NodeStoreSol<Double>::GetAlgSysDoublePointer()
  {
    ENTER_FCN("NodeStoreSol::GetAlgSysDoublePointer");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",
                            __FILE__,__LINE__);
#endif
    return data_.data_;
  }


  template<class TYPE>
  Double* NodeStoreSol<TYPE>::GetAlgSysDoublePointer()
  {
    ENTER_FCN("NodeStoreSol::GetAlgSysDoublePointer");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

    Error("Not implemented here",__FILE__,__LINE__);
    return NULL;
  }



  ///////// Transformation Operations ///////// 

  template<class TYPE>
  void NodeStoreSol<TYPE>::GetElemSolution(CFSVector & elemSol,
                                           const EntityIterator& it) const
  {
    ENTER_FCN( "NodeStoreSol::GetElemSolution" );

    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(elemSol);

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *result_, it );
    temp.Resize(eqns.GetSize());
    temp.Init();

     for( UInt i = 0; i < eqns.GetSize(); i++ ) {
       if ( eqns[i] != 0 ) {
         temp[i] = data_[(abs(eqns[i])-1)];
       } else {
         temp[i] = TYPE();
       }
     }
  }


  template<class TYPE>
  void NodeStoreSol<TYPE>::GetElemSolutionAsMatrix(CFSMatrix & elemSol, 
                                                   const EntityIterator& it) const
  {
    ENTER_FCN("NodeStoreSol::GetElemSolutionAsMatrix");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
    Matrix<TYPE> & temp = dynamic_cast<Matrix<TYPE>&>(elemSol);
    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *result_, it );

    UInt numFcns = (UInt) eqns.GetSize() / totalDofs_;    
    temp.Resize(totalDofs_, numFcns);

    for ( UInt iDof = 0; iDof < totalDofs_; iDof++ )
      for ( UInt iFcn=0; iFcn < numFcns; iFcn++ ) {
        Integer actEqn = eqns[iFcn * totalDofs_ + iDof];
        if ( actEqn != 0){
          temp[iDof][iFcn] = data_[(abs(actEqn)-1)];
        } else {
          temp[iDof][iFcn] = TYPE();
        }
      }
  }

  
  template<class TYPE>
  void NodeStoreSol<TYPE>::NodeSolutionToCoupling(CFSVector & couplingSol,
                                                  const StdVector<UInt>& nodeNumbers) const
  {
    ENTER_FCN("NodeStoreSol:::NodeSolutionToCoupling");
#ifdef CHECK_INITIALIZED
    if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(couplingSol);
    Integer eqnNr;

    temp.Resize(nodeNumbers.GetSize() * totalDofs_);
    for (UInt iNode=0; iNode<nodeNumbers.GetSize(); iNode++)
      for (UInt iDof=0; iDof<totalDofs_; iDof++)
        {
          eqnNr = eqnMap_->GetNodeEqn( nodeNumbers[iNode], iDof+1 );
          if (eqnNr != 0)
            temp.data_[iNode*totalDofs_ + iDof] = data_[abs((eqnNr)-1)];
          else
            temp.data_[iNode*totalDofs_ + iDof] = TYPE();
        }
  }

  
  //////////////////// Operators //////////////////
  template<class TYPE>
  NodeStoreSol<TYPE> & NodeStoreSol<TYPE>::operator= (const NodeStoreSol & x)
  {
    ENTER_FCN("NodeStoreSol::operator=(const NodeStoreSol &");
    Error("Not implemented here", __FILE__,__LINE__); 
    return *this;
  }

  template<class TYPE>
  BaseNodeStoreSol & NodeStoreSol<TYPE>::operator= (const BaseNodeStoreSol & x)
  {
    ENTER_FCN("NodeStoreSol::operator=(const BaseNodeStoreSol &");
    if ( &x == dynamic_cast<BaseNodeStoreSol*>(this))
      return (*this);

    const NodeStoreSol<TYPE> & temp = dynamic_cast<const NodeStoreSol<TYPE>&>(x);
  
    this->numNodes_ = temp.numNodes_;
    this->numSolutions_ = temp.numSolutions_;
    this->length_ = temp.length_;
    this->solTypes_ = temp.solTypes_;
    this->solOffset_ = temp.solOffset_;
    this->solDofs_ = temp.solDofs_; this->totalDofs_ = temp.totalDofs_;
    this->data_ = temp.data_;
    this->convertedData_ = temp.convertedData_;
 
    return dynamic_cast<BaseNodeStoreSol &> (*this);
  }

  template <class TYPE>
  void NodeStoreSol<TYPE>::Print( std::ostream& str ) {
    ENTER_IFCN( "NodeStoreSol::Print" );
    Error( "NodeStoreSol::Print() not implemented", __FILE__, __LINE__ );
  }

  // explicit template instantiation for SGI compiler
#ifdef __GNUC__
  template class NodeStoreSol<Double>;
  template class NodeStoreSol<Complex>;
#endif 

  // explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate NodeStoreSol<Double>
#pragma instantiate NodeStoreSol<Complex>
#endif

} //namespace
