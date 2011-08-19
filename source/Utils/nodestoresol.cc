// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "nodestoresol.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"

namespace CoupledField {


  // ***********************
  //   Default Constructor
  // ***********************
  template<class TYPE>
  NodeStoreSol<TYPE>::NodeStoreSol() {

  
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
    EXCEPTION( "Not implemented here" );
  }


  // ***************************
  //   Constructor (Version 2)
  // ***************************
  template<class TYPE>
  NodeStoreSol<TYPE>::NodeStoreSol( const UInt numNodes,
                                    const SolutionType solType,
                                    const UInt numDofs ) {
    EXCEPTION( "Not implemented here" );
  }


  // ***************************
  //   Constructor (Version 3)
  // ***************************
  template<class TYPE>
  NodeStoreSol<TYPE>::NodeStoreSol( const NodeStoreSol &x ) {
    ptGrid_ = x.ptGrid_;
    results_ = x.results_;
    elems_ = x.elems_;
    numNodes_ = x.numNodes_;
    numSolutions_ = x.numSolutions_;
    length_ = x.length_;
    lengthVector_ = x.lengthVector_;
    solTypes_ = x.solTypes_;
    solOffset_ = x.solOffset_;
    eqnOffset_ = x.eqnOffset_;
    solDofs_ = x.solDofs_;
    totalDofs_ = x.totalDofs_;
    eqnDofs_ = x.eqnDofs_;
    convertedData_ = x.convertedData_;

    data_ = x.data_;
  }


  // **************
  //   Destructor
  // **************
  template<class TYPE>
  NodeStoreSol<TYPE>::~NodeStoreSol() {
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
    EXCEPTION( "NodeStoreSol::Clear() not implemented" );
  }


  // ********
  //   Init
  // ********
  template<class TYPE>
  void NodeStoreSol<TYPE>::Init() {
    Init( TYPE() );
  }


  // *****************
  //   SetPtrEQNData
  // *****************
  template<class TYPE>
  void NodeStoreSol<TYPE>::SetPtrEQNData( shared_ptr<FeSpace> space,
                                          Grid *ptGrid ) {
    ptGrid_ = ptGrid;
    feSpace_ = space;
  }


  // ********
  //   Init
  // ********
  template<class TYPE>
  void NodeStoreSol<TYPE>::Init( const TYPE val ) {

  
#ifdef CHECK_INITIALIZED
    if (numSolutions_ == 0 || numNodes_ == 0 || solTypes_.size() == 0 ||
        solDofs_.size() == 0 ) {
      EXCEPTION( "NodeStoreSol::Init:"
                 << "\n numSolutions ....... " << numSolutions_
                 << "\n numNodes ........... " << numNodes_
                 << "\n totalDofs .......... " << totalDofs_
                 << "\n length ............. " << length_
                 << "\n size of solDofs .... " << solDofs_.size()
                 << "\n size of solTypes ... " << solTypes_.size()
                 << "\n size of offsets .... " << solOffset_.size()
                 << "\n Before calling Init(), the number of solutions, "
                 << "nodes, types and dofs has to be set to NON-ZERO value!" );
    }
#endif
  
    // only the first time the struct gets initialized
    //if ( length_ == 0 ) {
      if ( solDofs_.size() != numSolutions_ ||
           solTypes_.size() != numSolutions_ ) {
        EXCEPTION( "Inconsistent definition of Storesolution class. "
                   << "Maybe you have to call 'Clear()' before using a "
                   << "modified data layout!" );
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
      if ( feSpace_ == NULL ) {
        EXCEPTION("NodeStoreSol::Init: Pointer to EQN is NULL!" );
      }
      
      // Determine the length of the solution vector
      lengthVector_ = feSpace_->GetNumEquations();
      eqnDofs_ = 1;

      data_.Clear();
      data_.Resize(lengthVector_);
      data_.Init(val);
    //}
    // Assuming that the struct is already initialized
    // and it is only needed a change of the initial value
//    else{
//    	data_.Init(val);
//    }
    	
  }


  // *******************
  //   SetNumSolutions
  // *******************
  template<class TYPE>
  void NodeStoreSol<TYPE>::SetNumSolutions( const UInt nSols ) {
    numSolutions_ = nSols;
    length_ = 0;
  }


  // ***************
  //   SetNumNodes
  // ***************
  template<class TYPE>
  void NodeStoreSol<TYPE>::SetNumNodes( const UInt nNodes ) {
    numNodes_ = nNodes;
    length_ = 0;
  }


  // *******************
  //   SetSolutionType
  // *******************
  template<class TYPE>
  void NodeStoreSol<TYPE>::SetSolutionType( const SolutionType solType,
                                            const UInt numSolution ) {


#ifdef CHECK_INDEX
    if (numSolution >= numSolutions_) {
      EXCEPTION("NodeStoreSol: Index out of Bounds");
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
  
  

  template<class TYPE>
  void NodeStoreSol<TYPE>::SetSolutionTypeName( const SolutionType solType,
                                                const SolutionType solTypeName,
                                                const UInt numSolution ) {


#ifdef CHECK_INDEX
    if (numSolution >= numSolutions_) {
      EXCEPTION("NodeStoreSol: Index out of Bounds");
    }
#endif

    // Check, if the object contains only one entry and
    // if this only entry is overwritten
    if ( solTypes_.size() == 1 && numSolutions_ == 1 && numSolution == 0 ) {
      solTypes_.clear();
    }

    solTypes_[solType] = numSolution;
    solTypesName_[solTypeName] = numSolution;
    length_ = 0;
  }



  // **************
  //   SetNumDofs
  // **************
  template<class TYPE>
  void NodeStoreSol<TYPE>::SetNumDofs( const UInt dof,
                                       const SolutionType sol ) {


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


    solTypes.Resize(numSolutions_);
    std::map<SolutionType,UInt>::const_iterator it;
 
    for ( it = solTypes_.begin(); it != solTypes_.end(); it++ ) {
      solTypes[(*it).second] = (*it).first;
    }
  }

  template<class TYPE>
  TYPE& NodeStoreSol<TYPE>::operator()(UInt node, UInt dof)
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif

#ifdef CHECK_INDEX
    if (numSolutions_ > 1) 
      EXCEPTION("NodeStoreSol:operator(): Only defined objects with one type of solution!" );
#endif
    Integer eqnNr;
    eqnNr = feSpace_->GetNodeEqn( node+1, dof+1 );

    if (eqnNr > 0)
      return data_[abs(eqnNr)-1];
    else
      EXCEPTION("NodeStoreSol::operator(): This operator gives only writing access to non-BC nodes." );
  }

  template<class TYPE>
  TYPE  NodeStoreSol<TYPE>::operator()(UInt node, UInt dof) const
  {
  
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif

#ifdef CHECK_INDEX
    if (numSolutions_ > 1) 
      EXCEPTION("NodeStoreSol:operator(): Only defined objects with one type of solution!" );
#endif
  
    Integer eqnNr;
    eqnNr = feSpace_->GetNodeEqn( node+1,dof+1 );

    WARN ("Is this operator ever be used?");
    if (eqnNr > 0)
      return data_[abs(eqnNr)-1];
    else
      return TYPE();
  }


  template<class TYPE>
  void NodeStoreSol<TYPE>::GetGlobalSolVector(const SolutionType type,
                                              SingleVector & val) const
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) 
      EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif
    
    UInt globNumNodes = ptGrid_->GetNumNodes();
    StdVector<Integer> eqns;
    

    // Find related resultDof
    shared_ptr<ResultInfo> actResult;
    bool solFound = false;
    for( UInt i = 0; i < results_.GetSize(); i++ ) {
      if( results_[i]->resultType == type ) {
        actResult = results_[i];
        solFound = true;
        break;
      }
    }
    if( !solFound) {
      EXCEPTION("Solution type '" << SolutionTypeEnum.ToString(type)
                << "' is not contained in this solution object!");
    }

    UInt numDofs = actResult->dofNames.GetSize();
    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
    temp.Resize(globNumNodes *numDofs);
    temp.Init();

    // Iterate over all entities
    for( UInt i = 0; i < elems_.GetSize(); i++ ) {

      // Get iterator
      EntityIterator it = elems_[i]->GetIterator();

      for (it.Begin(); !it.IsEnd(); it++ ) {

        // Get connectivity
        StdVector<UInt> const & nodes  = it.GetElem()->connect;

        // Get related equation numbers
        //  ->GetEqns( eqns, *result_, it );

        // Iterate over all nodes
        for( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
          Integer globNode = nodes[iNode];

          // Iterate over all dofs
          for (UInt iDof = 0; iDof < numDofs; iDof++ ) {
            
            //WARNING!!: TODO This will cause a crash down in case of coupled PDEs
            Integer eqnNr = feSpace_->GetNodeEqn( globNode, iDof+1 );

            if( eqnNr != 0 ) {
              temp[(globNode-1)*numDofs + iDof] = data_[std::abs(eqnNr)-1];
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
  void NodeStoreSol<TYPE>::SetNodalResult(const UInt nodeNr, const SingleVector &val)
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif
  
#ifdef CHECK_INDEX
    if (nodeNr >= numNodes_) 
      EXCEPTION("NodeStoreSol::SetNodalResult(): index out of bounds" );
    if (val.GetSize() != totalDofs_)
      EXCEPTION("NodeStoreSol::SetNodalResult(): vector of incompatible dimension" );
#endif

    const Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
    for (UInt i=0; i<temp.GetSize(); i++) {
      UInt eqnNr = feSpace_->GetNodeEqn( nodeNr, i+1 );
      data_[eqnNr-1] = temp[i];
    }
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::GetNodalResult(const UInt nodeNr, SingleVector &val) const
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif
    EXCEPTION("Not implemented here" );
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::GetGlobalSolVectorSingleDof(const SolutionType type, 
                                                       const UInt dof, 
                                                       SingleVector & val) const
  {
#ifdef CHECK_INITIALIZED
    UInt solDof = (*solDofs_.find(type)).second;

    if (length_ == 0) 
      EXCEPTION("NodeStoreSol: Use of uninitialized object!" );


    std::string errMsg;
    if (dof > solDof)
      {
        EXCEPTION( "NodeStoreSol::GetSolVectorSingleDof: Desired dof "
                   << dof << " is higher than dofs of solution = "
                   << lexical_cast<std::string>(solDof) );
      }
#endif
  
    UInt offset = (*solOffset_.find(type)).second;
    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
    
    Integer eqnNr;
    //UInt globNumNodes = feSpace_->GetNumUnknowns();
    // TODO works only for Single PDE
    UInt numLocNodes = feSpace_->GetNumUnknowns();
    temp.Resize(numLocNodes);
    temp.Init();
    for (UInt iNode=0; iNode<numLocNodes; iNode++) {
      eqnNr = feSpace_->GetNodeEqn( iNode, dof+1+offset ); 
 
      if (eqnNr != 0) 
        temp[iNode-1] = data_[std::abs(eqnNr)-1];
      else if (eqnNr == 0)
        temp[iNode-1] = TYPE();
      else
        EXCEPTION( "Constraints not yet implemented!" );
    } 
  
  }  


  template<class TYPE>
  void NodeStoreSol<TYPE>::GetGlobalSolVectorSingleDof(const UInt dof, 
                                                       SingleVector & val) const{
#ifdef CHECK_INITIALIZED
    if (length_ == 0) {
      EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
    }
#endif
  
    //Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
 
    EXCEPTION ("Not working properly at the moment" );

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


#ifdef CHECK_INITIALIZED
    if (length_ == 0) {
      EXCEPTION( "NodeStoreSol::s): Use of uninitialized object!" );
    }
#endif

#ifdef CHECK_INDEX
    if ( numSolutions_ > 1 ) {
      EXCEPTION( "NodeStoreSol::Get(): Only used for single solution objects!" );
    }
#endif

    Integer eqnNr;
    eqnNr = feSpace_->GetNodeEqn( nodeNr + 1, dof + 1 );
    if ( eqnNr != 0 ) {
      ret = data_[ abs(eqnNr)- 1];
    }
    else {
      ret =  TYPE(0.0);
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
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!");
#endif

    // Find related resultDof
    shared_ptr<ResultInfo> actResult;
    bool solFound = false;
    for( UInt i = 0; i < results_.GetSize(); i++ ) {
      if( results_[i]->resultType == type ) {
        actResult = results_[i];
        solFound = true;
        break;
      }
    }
    if( !solFound) {
      EXCEPTION("Solution type '" << SolutionTypeEnum.ToString(type)
                << "' is not contained in this solution object!");
    }  

    Integer eqnNr;
    eqnNr = feSpace_->GetNodeEqn( nodeNr+1, dof+1);

    if (eqnNr != 0)
      ret = data_[abs(eqnNr)-1];
    else
      ret =  TYPE();
  }


  template<class TYPE>
  void NodeStoreSol<TYPE>::Set(const SolutionType type, const UInt nodeNr, const UInt dof, const TYPE val)
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif

    UInt offset = (*solOffset_.find(type)).second;

    Integer eqnNr;
  
    eqnNr = feSpace_->GetNodeEqn( nodeNr+1, offset+dof+1 );

    if (eqnNr != 0)
      data_[abs(eqnNr)-1] = val;
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::Set(const Integer eqnNr, 
                               const TYPE val)
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!");
#endif

    if (eqnNr != 0)
      data_[abs(eqnNr-1)] = val;
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::Add(const SolutionType type, const UInt nodeNr, const UInt dof, const TYPE val) const
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif
    EXCEPTION("Not implemented here" );
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::Add(const Integer eqnNr, 
                               const TYPE val)
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!");
#endif
    if (eqnNr != 0)
      data_[abs(eqnNr-1)] += val;
  }


  template<class TYPE>
  void NodeStoreSol<TYPE>::SetAlgSysVector(const SingleVector & val)
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) WARN ("Use of uninitialized object!");
#endif
#ifdef CHECK_INDEX
   //  if (val.GetSize() !=  lengthVector_)
//       {
//         std::cerr << "Vector has Size" << lengthVector_ << std::endl;
//         std::cerr << "Your vector has size " << val.GetSize() << std::endl;
//         EXCEPTION("NodeStoreSol::SetAlgSysVector(): Vector has wrong size!" );
//       }
#endif

    const  Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
    data_ = temp;
  }
  
  template<class TYPE> 
  void NodeStoreSol<TYPE>::GetAlgSysVector(SingleVector & val) const
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif
  
    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
    temp = data_;
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::GetAlgSysVectorPointer(SingleVector* &ptrToVec)
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif

    ptrToVec = (SingleVector*) &data_;
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::GetAlgSysVectorPointer(Vector<TYPE>* &ptrToVec)
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif

    ptrToVec =  &data_;
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::CopyFromAlgSysDataPointer(Double * ptr)
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
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
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif
    for (UInt i=0; i<lengthVector_; i++)
      {
        data_[i] = Complex(ptr[2*i],ptr[2*i+1]);
      }

  }


  template<class TYPE>
  void NodeStoreSol<TYPE>::SetAlgSysDataPointer( UInt size, TYPE * ptr)
  {
  
    data_.Replace( size, ptr, false );

  }


  template<>
  Double* NodeStoreSol<Double>::GetAlgSysDoublePointer()
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif
    return data_.data_;
  }


  template<class TYPE>
  Double* NodeStoreSol<TYPE>::GetAlgSysDoublePointer()
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );
#endif

    EXCEPTION("Not implemented here");
    return NULL;
  }



  ///////// Transformation Operations ///////// 

  template<class TYPE>
  void NodeStoreSol<TYPE>::GetElemSolution(SingleVector & elemSol,
                                           const EntityIterator& it, 
                                           UInt solIndex ) const
  {

    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(elemSol);

    StdVector<Integer> eqns;
    //eqnMap_->GetEqns( eqns, *results_[solIndex], it );
    feSpace_->GetEqns(eqns, it);
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
  void NodeStoreSol<TYPE>::GetElemSolutionAsMatrix(DenseMatrix & elemSol, 
                                                   const EntityIterator& it,
                                                   UInt solIndex, const SingleVector* ext_data) const
  {
#ifdef CHECK_INITIALIZED
    if(ext_data != NULL) 
    { 
      if(ext_data->GetSize() == 0) EXCEPTION("External data uninitialized"); 
      if(ext_data->GetSize() != data_.GetSize()) EXCEPTION("External data does not match internal one"); 
    } 
    else 
    {
      if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!" );      
    }
#endif
  
    // choose either external data or if not given use own one 
    const Vector<TYPE>& data = ext_data == NULL ? data_ : dynamic_cast<const Vector<TYPE>&>(*ext_data);  
    
    Matrix<TYPE> & temp = dynamic_cast<Matrix<TYPE>&>(elemSol);
    StdVector<Integer> eqns;
    //eqnMap_->GetEqns( eqns, *results_[solIndex], it );

    UInt numDofs =  results_[solIndex]->dofNames.GetSize();
    UInt numFcns = (UInt) eqns.GetSize() / numDofs;
    temp.Resize( numDofs, numFcns);

    for ( UInt iDof = 0; iDof < numDofs; iDof++ )
      for ( UInt iFcn=0; iFcn < numFcns; iFcn++ ) {
        Integer actEqn = eqns[iFcn * numDofs + iDof];
        if ( actEqn != 0){
          temp[iDof][iFcn] = data[(abs(actEqn)-1)];
        } else {
          temp[iDof][iFcn] = TYPE();
        }
      }
  }

  template<class TYPE>
  void NodeStoreSol<TYPE>::NodeSolutionToCoupling(SingleVector & couplingSol,
                                                  const StdVector<UInt>& nodeNumbers,
                                                  UInt solIndex) const
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("NodeStoreSol: Use of uninitialized object!");
#endif

    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(couplingSol);
    Integer eqnNr;

    UInt numDofs =  results_[solIndex]->dofNames.GetSize();
    temp.Resize(nodeNumbers.GetSize() * numDofs);
    for (UInt iNode=0; iNode<nodeNumbers.GetSize(); iNode++)
      for (UInt iDof=0; iDof<numDofs; iDof++)
        {
          eqnNr = feSpace_->GetNodeEqn( nodeNumbers[iNode], iDof+1 );
          if (eqnNr != 0)
            temp.data_[iNode*totalDofs_ + iDof] = data_[abs(eqnNr)-1];
          else
            temp.data_[iNode*numDofs + iDof] = TYPE();
        }
  }


  //////////////////// Operators //////////////////
  template<class TYPE>
  NodeStoreSol<TYPE> & NodeStoreSol<TYPE>::operator= (const NodeStoreSol & x)
  {
    EXCEPTION("Not implemented here" );
    return *this;
  }

  template<class TYPE>
  BaseNodeStoreSol & NodeStoreSol<TYPE>::operator= (const BaseNodeStoreSol & x)
  {
    if ( &x == dynamic_cast<BaseNodeStoreSol*>(this))
      return (*this);

    const NodeStoreSol<TYPE> & temp = dynamic_cast<const NodeStoreSol<TYPE>&>(x);
  
    this->data_ = temp.data_;

    this->ptGrid_ = temp.ptGrid_;
    this->results_ = temp.results_;
    this->elems_ = temp.elems_;
    this->numNodes_ = temp.numNodes_;
    this->numSolutions_ = temp.numSolutions_;
    this->length_ = temp.length_;
    this->lengthVector_ = temp.lengthVector_;
    this->solTypes_ = temp.solTypes_;
    this->solTypesName_ = temp.solTypesName_;
    this->solOffset_ = temp.solOffset_;
    this->solDofs_ = temp.solDofs_;
    this->totalDofs_ = temp.totalDofs_;
    this->eqnDofs_ = temp.eqnDofs_;
    this->convertedData_ = temp.convertedData_;
 
    return dynamic_cast<BaseNodeStoreSol &> (*this);
  }

  template <class TYPE>
  void NodeStoreSol<TYPE>::Print( std::ostream& str ) {
    EXCEPTION( "NodeStoreSol::Print() not implemented" );
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
