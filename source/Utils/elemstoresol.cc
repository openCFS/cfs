// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "elemstoresol.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField{



  template<class TYPE>
  ElemStoreSol<TYPE>::ElemStoreSol()
  {
  
    numElems_ = 0;
    numSolutions_ = 0;
  }
  
  template<class TYPE>
  ElemStoreSol<TYPE>::ElemStoreSol(UInt numElems, 
                                   StdVector<SolutionType> solTypes, 
                                   StdVector<UInt> solDofs)
  {
    EXCEPTION("Not implemented here");
  }

  template<class TYPE>
  ElemStoreSol<TYPE>::ElemStoreSol(const UInt numElems,
                                   const SolutionType solType,
                                   const UInt numDofs)
  {
    EXCEPTION("Not implemented here");
  }


  template<class TYPE>
  ElemStoreSol<TYPE>::ElemStoreSol(const ElemStoreSol & x) 
  {
    EXCEPTION("Not implemented here");
  }



  template<class TYPE>
  ElemStoreSol<TYPE>::~ElemStoreSol() 
  {
 
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::SetPtrEQNData( EqnMap * eqnMap,
                                          Grid * ptGrid)
  {
    eqnMap_ = eqnMap;
    ptGrid_ = ptGrid;
  }


  template<class TYPE>
  void ElemStoreSol<TYPE>::Clear()
  {
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::Init()  
  {
    Init(TYPE());
  }


  template<class TYPE>
  void ElemStoreSol<TYPE>::Init(const TYPE val)  
  {

#ifdef CHECK_INITIALIZED
    if (numSolutions_ == 0 || numElems_ == 0 \
        || solTypes_.size() == 0 || solDofs_.size() == 0) 
      {
        std::cerr << "Error in StoreSl::Init():" << std::endl;
        std::cerr << "NumSolutions: " << numSolutions_ << std::endl;
        std::cerr << "NumElems: " << numElems_ << std::endl;
        std::cerr << "TotalDofs: " << totalDofs_ << std::endl;
        std::cerr << "Length: " << length_ << std::endl;
        std::cerr << "Size of solDofs: " << solDofs_.size() << std::endl;
        std::cerr << "Size of solTypes: " << solTypes_.size() << std::endl;
        std::cerr << "Size of offsets: " << solOffset_.size() << std::endl;
        EXCEPTION("ElemStoreSol::Init(): Before calling Init(), the number of solutions,\
           elems, types and dofs has to be set to NONZERO value!");
      }
  
#endif
  
    // only the first time the struct gets initialized
    if (length_ == 0)
      {
 
        if (solDofs_.size() != numSolutions_ || solTypes_.size() != numSolutions_)
          EXCEPTION("Inconsistent definition of Storesolution class.\
                     Eventually you have to call 'Clear()' before using a modified data layout!");
      
        totalDofs_ = 0;
      
        // iterate over all solutiontypes and
        // sum up their number of dof
        std::map<SolutionType,UInt>::iterator it;
        for (it = solDofs_.begin(); it!=solDofs_.end(); it++)
          {
            // set offset of current solution
            // w.r.t. to starting position
            solOffset_[(*it).first] = totalDofs_;
            totalDofs_ += (*it).second;     
          }
      
        length_ = totalDofs_ * numElems_;
        data_.Resize(length_);
      }

    data_.Init(val);
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::SetNumSolutions(const UInt nSols)
  {

    numSolutions_ = nSols;
    length_ = 0;
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::SetNumElems(const UInt nElems)
  {
    numElems_ = nElems;
    length_ = 0;
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::SetSolutionType(const SolutionType solType, const UInt numSolution)
  {

#ifdef CHECK_INDEX
    if (numSolution >= numSolutions_) EXCEPTION("ElemStoreSol: Index out of Bounds");
#endif

    solTypes_[solType] = numSolution;
    length_ = 0;
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::SetNumDofs(const UInt dof, const SolutionType sol)
  {
  
    // check if only one dof was assigned
    // -> only one entry exists
    if (numSolutions_ == 1 && sol == NO_SOLUTION_TYPE)
      solDofs_[(*solTypes_.begin()).first] = dof;
    else
      solDofs_[sol] = dof;

    length_ = 0;
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::GetSolutionTypes(StdVector<SolutionType> &solTypes) const
  {
  
    solTypes.Resize(numSolutions_);
    std::map<SolutionType,UInt>::const_iterator it;
 
    for (it = solTypes_.begin(); it!=solTypes_.end(); it++)
      solTypes[(*it).second] = (*it).first;
  
  }

  template<class TYPE>
  TYPE& ElemStoreSol<TYPE>::operator()(UInt node, UInt dof)
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("ElemStoreSol: Use of uninitialized object!");
#endif

#ifdef CHECK_INDEX
    if (numSolutions_ > 1) 
      EXCEPTION("ElemStoreSol:operator(): Only defined objects with one type of solution!");
#endif

    return data_[node * totalDofs_ + dof];
  }

  template<class TYPE>
  TYPE  ElemStoreSol<TYPE>::operator()(UInt node, UInt dof) const
  {
  
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("ElemStoreSol: Use of uninitialized object!");
#endif

#ifdef CHECK_INDEX
    if (numSolutions_ > 1) 
      EXCEPTION("ElemStoreSol:operator(): Only defined objects with one type of solution!");
#endif
 
    return data_[node * totalDofs_ + dof];
  }


  template<class TYPE>
  void ElemStoreSol<TYPE>::SetElemResult(const UInt elemNr, const SingleVector &val)
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("ElemStoreSol: Use of uninitialized object!");
#endif
  
#ifdef CHECK_INDEX
    if (elemNr >= numElems_) 
      EXCEPTION("ElemStoreSol::SetElemResult(): index out of bounds");
    if (val.GetSize() != totalDofs_)
      EXCEPTION("ElemStoreSol::SetElemResult(): vector of incompatible dimension");
#endif

    const Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
    for (UInt i=0; i<temp.GetSize(); i++){
      data_[elemNr*totalDofs_ + i] = temp[i];
      // std::cout<<"data_"<< data_[elemNr*totalDofs_ + i] << std::endl;
    }
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::GetElemResult(const UInt elemNr, SingleVector &val) const
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("ElemStoreSol: Use of uninitialized object!");
#endif
    
   Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
   for (UInt i=0; i < totalDofs_; i++){
      temp[i] = data_[elemNr*totalDofs_ + i];
      // std::cout<<"data_"<< data_[elemNr*totalDofs_ + i] << std::endl;
    }
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::GetGlobalSolVector(const SolutionType solType, SingleVector & val) const
  {
    // killme ! the solType is not queried
#ifdef CHECK_INITIALIZED
    if (length_ == 0) WARN("Use of uninitialized object!");
#endif
  
    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
    temp.Resize(ptGrid_->GetNumElems()*totalDofs_);
    temp.Init();

    // Loop over all PDE elements
    for (UInt iElem=1; iElem<numElems_+1; iElem++)
      // Loop over all dimensions
      for (UInt iDof=0; iDof<totalDofs_; iDof++)
        {
          //temp.data_[(mapping_[iElem]-1)*totalDofs_ + iDof] = data_[iElem*totalDofs_ + iDof];
          // temp.data_[(ptEQN_->PDE2MeshElem(iElem)-1)*totalDofs_ + iDof] = 
//             data_[(iElem-1)*totalDofs_ + iDof];
          temp.data_[(eqnMap_->Pde2MeshElem(iElem)-1)*totalDofs_ + iDof] = 
            data_[(iElem-1)*totalDofs_ + iDof];
        }


  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::GetGlobalSolVectorSingleDof(const SolutionType type, const UInt dof, SingleVector & val) const
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("ElemStoreSol: Use of uninitialized object!");
#endif
    EXCEPTION("Not implemented here");
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::GetGlobalSolVectorSingleDof(const UInt dof, SingleVector & val) const
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("ElemStoreSol: Use of uninitialized object!");
#endif

    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
    temp.Resize( ptGrid_->GetNumElems() );
    temp.Init();

    // Loop over all PDE elements
    for (UInt iElem=1; iElem<numElems_+1; iElem++)
      // Loop over all dimensions
      temp.data_[eqnMap_->Pde2MeshElem(iElem)-1] = data_[(iElem-1)*totalDofs_ + dof];
  
  }


  template<class TYPE>
  void ElemStoreSol<TYPE>::Get(const UInt elemNr, const UInt dof, TYPE & ret) const
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("ElemStoreSol: Use of uninitialized object!");
#endif

#ifdef CHECK_INDEX
    if (numSolutions_ > 1)
      EXCEPTION("ElemStoreSol::Get(): Only used for single solution objects!");
    if (elemNr > numElems_)
      EXCEPTION("ElemStoreSol::Get(): Index out of bounds");
#endif

    ret = data_[elemNr * totalDofs_ + dof]; 
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::Get(const SolutionType type, const UInt elemNr, const UInt dof, TYPE & ret) const
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("ElemStoreSol: Use of uninitialized object!");
#endif
    UInt offset = (*solOffset_.find(type)).second;

    ret = data_[elemNr * totalDofs_ + offset + dof];
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::Set(const SolutionType type, const UInt elemNr, const UInt dof, const TYPE val)
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("ElemStoreSol: Use of uninitialized object!");
#endif
    EXCEPTION("Not implemented here");
  }

  template<class TYPE>
  void ElemStoreSol<TYPE>::Add(const SolutionType type, const UInt elemNr, const UInt dof, const TYPE val) const
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("ElemStoreSol: Use of uninitialized object!");
#endif
    EXCEPTION("Not implemented here");
  }

  ///////// Transformation Operations ///////// 

  template<class TYPE>
  void ElemStoreSol<TYPE>::TransformElemSolution(SingleVector & transformedSolution,
                                                 Grid * ptGrid) const
  {
#ifdef CHECK_INITIALIZED
    if (length_ == 0) WARN("Use of uninitialized object!");
#endif
  
    Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(transformedSolution);
    temp.Resize(ptGrid->GetNumElems()*totalDofs_);
    temp.Init();

    // Loop over all PDE elements
    for (UInt iElem=1; iElem<numElems_+1; iElem++)
      // Loop over all dimensions
      for (UInt iDof=0; iDof<totalDofs_; iDof++)
        {
          //temp.data_[(mapping_[iElem]-1)*totalDofs_ + iDof] = data_[iElem*totalDofs_ + iDof];
          temp.data_[(eqnMap_->Pde2MeshElem(iElem)-1)*totalDofs_ + iDof] = 
            data_[(iElem-1)*totalDofs_ + iDof];
        }
  }


  template<class TYPE> 
  void ElemStoreSol<TYPE>::ElemSolutionToCoupling(SingleVector & couplingSol,
                                                  const StdVector<Elem*>& elements,
                                                  const SingleVector & elemSol) const
  {

#ifdef CHECK_INITIALIZED
    if (length_ == 0) EXCEPTION("ElemStoreSol: Use of uninitialized object!");
#endif

    EXCEPTION("Not implemented here"); 

  }
  
  
  //////////////////// Operators //////////////////
  template<class TYPE>
  ElemStoreSol<TYPE> & ElemStoreSol<TYPE>::operator= (const ElemStoreSol & x)
  {
    EXCEPTION("Not implemented here"); 
  }

  template<class TYPE>
  BaseElemStoreSol & ElemStoreSol<TYPE>::operator= (const BaseElemStoreSol & x)
  {
    if ( &x == dynamic_cast<ElemStoreSol*>(this))
      return (*this);

    const ElemStoreSol<TYPE> & temp = dynamic_cast<const ElemStoreSol<TYPE>&>(x);
  
    this->numElems_ = temp.numElems_;
    this->numSolutions_ = temp.numSolutions_;
    this->length_ = temp.length_;
    this->solTypes_ = temp.solTypes_;
    this->solOffset_ = temp.solOffset_;
    this->solDofs_ = temp.solDofs_;
    this->totalDofs_ = temp.totalDofs_;
    this->data_ = temp.data_;
    this->convertedData_ = temp.convertedData_;
 
    return dynamic_cast<BaseElemStoreSol &> (*this);
  }

  template <class TYPE>
  void ElemStoreSol<TYPE>::Print( std::ostream& str ) {
    EXCEPTION("ElemStoreSol::Print() not implemented");
  }

  // explicit template instantiation for GCC compiler
#ifdef __GNUC__
  template class ElemStoreSol<Double>;
  template class ElemStoreSol<Complex>;
#endif 


  // explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate ElemStoreSol<Double>
#pragma instantiate ElemStoreSol<Complex>
#endif

} //namespace
