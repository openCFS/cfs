#include "elemstoresol.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField{



template<class TYPE>
ElemStoreSol<TYPE>::ElemStoreSol()
{
  ENTER_FCN( "ElemStoreSol::ElemStoreSol()" );
  
  numNodes_ = 0;
  numSolutions_ = 0;
  length_ = 0;
  totalDofs_ = 0;

}
  
template<class TYPE>
ElemStoreSol<TYPE>::ElemStoreSol(Integer numNodes, 
			 StdVector<SolutionType> solTypes, 
			 StdVector<Integer> solDofs)
{
  ENTER_FCN( "ElemStoreSol::ElemStoreSol(numNodes, solTypes, solDofs" );
  Error( "Not implemented here", __FILE__, __LINE__ );
}

template<class TYPE>
ElemStoreSol<TYPE>::ElemStoreSol(const Integer numNodes,
			 const SolutionType solType,
			 const Integer numDofs)
{
  ENTER_FCN( "ElemStoreSol::ElemStoreSol(numNodes, solType, soDofs" );
  Error( "Not implemented here", __FILE__, __LINE__ );
}


template<class TYPE>
ElemStoreSol<TYPE>::ElemStoreSol(const ElemStoreSol & x) 
{
  ENTER_FCN( "ElemStoreSol::ElemStoreSol(const ElemStoreSol)" );
  Error( "Not implemented here", __FILE__, __LINE__ );
}



template<class TYPE>
ElemStoreSol<TYPE>::~ElemStoreSol() 
{
  ENTER_FCN( "ElemStoreSol::~ElemStoreSol" );
 
  
}

template<class TYPE>
void ElemStoreSol<TYPE>::Clear()
{
  ENTER_FCN( "ElemStoreSol::Clear()" );
}

template<class TYPE>
void ElemStoreSol<TYPE>::Init()  
{
  Init(TYPE());
}


template<class TYPE>
void ElemStoreSol<TYPE>::Init(const TYPE val)  
{
  ENTER_FCN( "ElemStoreSol::Init" );

#ifdef CHECK_INITIALIZED
  if (numSolutions_ == 0 || numNodes_ == 0 \
     || solTypes_.size() == 0 || solDofs_.size() == 0) 
    {
      std::cerr << "Error in StoreSl::Init():" << std::endl;
      std::cerr << "NumSolutions: " << numSolutions_ << std::endl;
      std::cerr << "NumNodes: " << numNodes_ << std::endl;
      std::cerr << "TotalDofs: " << totalDofs_ << std::endl;
      std::cerr << "Length: " << length_ << std::endl;
      std::cerr << "Size of solDofs: " << solDofs_.size() << std::endl;
      std::cerr << "Size of solTypes: " << solTypes_.size() << std::endl;
      std::cerr << "Size of offsets: " << solOffset_.size() << std::endl;
      Error("ElemStoreSol::Init(): Before calling Init(), the number of solutions,\
           nodes, types and dofs has to be set to NONZERO value!",__FILE__,__LINE__);
    }
  
#endif
  
  // only the first time the struct gets initialized
  if (length_ == 0)
    {
 
      if (solDofs_.size() != numSolutions_ || solTypes_.size() != numSolutions_)
	Error("Inconsistent definition of Storesolution class.\
                     Eventually you have to call 'Clear()' before using a modified data layout!",
		    __FILE__, __LINE__);
      
      totalDofs_ = 0;
      
      // iterate over all solutiontypes and
      // sum up their number of dof
      std::map<SolutionType,Integer>::iterator it;
      for (it = solDofs_.begin(); it!=solDofs_.end(); it++)
	{
      // set offset of current solution
      // w.r.t. to starting position
	  solOffset_[(*it).first] = totalDofs_;
	  totalDofs_ += (*it).second;	  
	}
      
      length_ = totalDofs_ * numNodes_;
      data_.Resize(length_);
    }

  data_.Init(val);
}

template<class TYPE>
void ElemStoreSol<TYPE>::SetNumSolutions(const Integer nSols)
{
  ENTER_IFCN("ElemStoreSol::SetNumSolutions");

  numSolutions_ = nSols;
  length_ = 0;
}

template<class TYPE>
void ElemStoreSol<TYPE>::SetNumNodes(const Integer nNodes)
{
  ENTER_IFCN("ElemStoreSol::SetNumNodes");
  numNodes_ = nNodes;
  length_ = 0;
}

template<class TYPE>
void ElemStoreSol<TYPE>::SetSolutionType(const SolutionType solType, const Integer numSolution)
{
  ENTER_IFCN("ElemStoreSol::SetSolutionType");

#ifdef CHECK_INDEX
  if (numSolution >= numSolutions_) Error("ElemStoreSol: Index out of Bounds",__FILE__,__LINE__);
#endif

  solTypes_[solType] = numSolution;
  length_ = 0;
}

template<class TYPE>
void ElemStoreSol<TYPE>::SetNumDofs(const Integer dof, const SolutionType sol)
{
  ENTER_IFCN("ElemStoreSol::SetDof");
  
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
  ENTER_FCN( "ElemStoreSol::GetSolutionTypes" );
  
  solTypes.Resize(numSolutions_);
  std::map<SolutionType,Integer>::const_iterator it;
 
    for (it = solTypes_.begin(); it!=solTypes_.end(); it++)
      solTypes[(*it).second] = (*it).first;
  
}

template<class TYPE>
TYPE& ElemStoreSol<TYPE>::operator()(Integer node, Integer dof)
{
  ENTER_IFCN("ElemStoreSol::operator()");
 #ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (numSolutions_ > 1) 
    Error("ElemStoreSol:operator(): Only defined objects with one type of solution!",__FILE__,__LINE__);
#endif

  return data_[node * totalDofs_ + dof];
}

template<class TYPE>
TYPE  ElemStoreSol<TYPE>:: operator()(Integer node, Integer dof) const
{
  ENTER_IFCN("ElemStoreSol::operator()");
  
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (numSolutions_ > 1) 
    Error("ElemStoreSol:operator(): Only defined objects with one type of solution!",__FILE__,__LINE__);
#endif
 
  return data_[node * totalDofs_ + dof];
}


template<class TYPE>
void ElemStoreSol<TYPE>::GetSolVector(const SolutionType type, CFSVector & val) const
{
  ENTER_FCN("ElemStoreSol::GetSolVector");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  Integer offset = (*solOffset_.find(type)).second;
  Integer dof = (*solDofs_.find(type)).second;
  val.Resize(dof*numNodes_);

  Vector<TYPE> & ret = dynamic_cast<Vector<TYPE>&>(val);
  for (Integer iNode=0; iNode<numNodes_; iNode++)
    for (Integer iDof=0; iDof<dof; iDof++)
      ret[iNode*dof+iDof] = data_[iNode*totalDofs_+iDof+offset];
}

template<class TYPE>
void ElemStoreSol<TYPE>::SetSolVector(const SolutionType type, const CFSVector & val)
{
  ENTER_FCN("ElemStoreSol::SetSolVector");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif 
  
#ifdef CHECK_INDEX
  if (val.GetSize() != data_.GetSize())
    Error("ElemStoreSol::SetSolVector(): Vector has incompatible dimensions!",__FILE__,__LINE__);
#endif

  const Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
  for (Integer i=0; i<temp.GetSize(); i++)
    data_[i] = temp[i];
}

template<class TYPE>
void ElemStoreSol<TYPE>::GetSolution(const SolutionType type, BaseStoreSol & val) const
{
  ENTER_FCN("ElemStoreSol::GetSolution");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Integer offset = (*solOffset_.find(type)).second;
  Integer dof = (*solDofs_.find(type)).second;
    
  ElemStoreSol<TYPE> & temp = dynamic_cast<ElemStoreSol<TYPE>&>(val);

  // delete old map
  temp.solDofs_.clear();
  temp.solTypes_.clear();
  temp.solOffset_.clear();

  temp.numNodes_ = numNodes_;
  temp.numSolutions_ = 1;
  temp.solTypes_[type] = 0;
  temp.solOffset_[type] = 0;
  temp.solDofs_[type] = dof;
  temp.totalDofs_ = dof;
  temp.length_ = temp.numNodes_ * totalDofs_;
  temp.data_.Resize(temp.length_);

  for (Integer iNode=0; iNode<numNodes_; iNode++)
    for (Integer iDof=0; iDof<dof; iDof++)
      temp.data_[iNode*dof+iDof] = data_[iNode*totalDofs_+iDof+offset];
}

template<class TYPE>
void ElemStoreSol<TYPE>::SetNodalResult(const Integer nodeNr, const CFSVector &val)
{
  ENTER_FCN("ElemStoreSol::SetNodalResult");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
#ifdef CHECK_INDEX
  if (nodeNr >= numNodes_) 
    Error("ElemStoreSol::SetNodalResult(): index out of bounds",__FILE__,__LINE__);
  if (val.GetSize() != totalDofs_)
    Error("ElemStoreSol::SetNodalResult(): vector of incompatible dimension",__FILE__,__LINE__);
#endif

  const Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
  for (Integer i=0; i<temp.GetSize(); i++)
    data_[nodeNr*totalDofs_ + i] = temp[i];
}

template<class TYPE>
void ElemStoreSol<TYPE>::GetNodalResult(const Integer nodeNr, CFSVector &val) const
{
  ENTER_FCN("ElemStoreSol::GetNodalResult");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Error("Not implemented here", __FILE__,__LINE__);
}

template<class TYPE>
void ElemStoreSol<TYPE>::GetSolVectorSingleDof(const SolutionType type, const Integer dof, CFSVector & val) const
{
  ENTER_FCN("ElemStoreSol::GetSolVectorSingleDof");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Error("Not implemented here", __FILE__,__LINE__);
}

template<class TYPE>
void ElemStoreSol<TYPE>::GetSolVectorSingleDof(const Integer dof, CFSVector & val) const
{
  ENTER_FCN("ElemStoreSol::GetSolVectorSingleDof");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Error("Not implemented here", __FILE__,__LINE__);
}


template<class TYPE>
void ElemStoreSol<TYPE>::Get(const Integer nodeNr, const Integer dof, TYPE & ret) const
{
  ENTER_FCN("ElemStoreSol::Get");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (numSolutions_ > 1)
    Error("ElemStoreSol::Get(): Only used for single solution objects!",__FILE__,__LINE__);
  if (nodeNr > numNodes_)
    Error("ElemStoreSol::Get(): Index out of bounds",__FILE__,__LINE__);
#endif

  ret = data_[nodeNr * totalDofs_ + dof]; 
}

template<class TYPE>
void ElemStoreSol<TYPE>::Get(const SolutionType type, const Integer nodeNr, const Integer dof, TYPE & ret) const
{
  ENTER_FCN("ElemStoreSol::Get");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Integer offset = (*solOffset_.find(type)).second;

  ret = data_[nodeNr * totalDofs_ + offset + dof];
}

template<class TYPE>
void ElemStoreSol<TYPE>::Set(const SolutionType type, const Integer nodeNr, const Integer dof, const TYPE val)
{
  ENTER_FCN("ElemStoreSol::Set");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Error("Not implemented here", __FILE__,__LINE__);
}

template<class TYPE>
void ElemStoreSol<TYPE>::Add(const SolutionType type, const Integer nodeNr, const Integer dof, const TYPE val) const
{
  ENTER_FCN("ElemStoreSol::Set");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Error("Not implemented here", __FILE__,__LINE__);
}

template<class TYPE>
void ElemStoreSol<TYPE>::SetCompleteVector(const CFSVector & val)
{
  ENTER_FCN("ElemStoreSol::SetCompleteVector");
#ifdef CHECK_INITIALIZED
   if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
#ifdef CHECK_INDEX
   if (val.GetSize() !=  length_)
     Error("ElemStoreSol::SetCompleteVector(): Vector has wrong size!",__FILE__,__LINE__);
#endif

   const  Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
   data_ = temp;
}
  
template<class TYPE> 
void ElemStoreSol<TYPE>::GetCompleteVector(CFSVector & val) const
{
 ENTER_FCN("ElemStoreSol::SetCompleteVector");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
  Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
  temp = data_;
}

template<class TYPE>
void ElemStoreSol<TYPE>::GetVectorPointer(CFSVector* &ptrToVec)
{
  ENTER_FCN("ElemStoreSol::GetVectorPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  ptrToVec = (CFSVector*) &data_;
}

template<class TYPE>
void ElemStoreSol<TYPE>::GetVectorPointer(Vector<TYPE>* &ptrToVec)
{
  ENTER_FCN("ElemStoreSol::GetVectorPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  ptrToVec =  &data_;
}

template<class TYPE>
void ElemStoreSol<TYPE>::CopyFromDataPointer(TYPE * ptr)
{
  ENTER_FCN("ElemStoreSol::CopyFromDataPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  for (Integer i=0; i<lengthVector_; i++)
    {
      data_[i] = ptr[i];
      //std::cerr << "Local Node [" << i+1 <<"] = " << ptr[i] << std::endl;
    }

}

template<class TYPE>
void ElemStoreSol<TYPE>::SetDataPointer(TYPE * ptr)
{
  ENTER_FCN("ElemStoreSol::SetDataPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",
				__FILE__,__LINE__);
#endif
  
  data_.data_ = ptr;

}

template<class TYPE>
void ElemStoreSol<TYPE>::GetDataPointer(TYPE* &ptr)
{
  ENTER_FCN("ElemStoreSol::GetDataPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",
				__FILE__,__LINE__);
#endif
 
  ptr =  data_.data_;
}


template<>
Double* ElemStoreSol<Double>::GetDoublePointer()
{
  ENTER_FCN("ElemStoreSol::GetDoublePointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",
				__FILE__,__LINE__);
#endif
  return data_.data_;
}


template<class TYPE>
Double* ElemStoreSol<TYPE>::GetDoublePointer()
{
  ENTER_FCN("ElemStoreSol::GetDoublePointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  Error("Not implemented here",__FILE__,__LINE__);
}



///////// Transformation Operations ///////// 

template<class TYPE>
void ElemStoreSol<TYPE>::TransformElemSolution(CFSVector & transformedSolution,
					   Grid * ptGrid,
					   const Integer level) const
{
  ENTER_FCN("ElemStoreSol::TransformElemSolution");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Warning("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
  Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(transformedSolution);
  temp.Resize(ptGrid->GetMaxnumElem(level)*totalDofs_);

  // Loop over all PDE elements
  for (Integer iElem=0; iElem<mapping_.GetSize(); iElem++)
    // Loop over all dimensions
    for (Integer iDof=0; iDof<totalDofs_; iDof++)
      {
      temp.data_[(mapping_[iElem]-1)*totalDofs_ + iDof] = data_[iElem*totalDofs_ + iDof];
      }
}


template<class TYPE> 
void ElemStoreSol<TYPE>::ElemSolutionToCoupling(CFSVector & couplingSol,
					    const StdVector<Elem*>& elements,
					    const CFSVector & elemSol) const
{
  ENTER_FCN("ElemStoreSol::ElemSolutionToCoupling");

#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("ElemStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  Error("Not implemented here", __FILE__,__LINE__); 

}
  
  
//////////////////// Operators //////////////////
template<class TYPE>
ElemStoreSol<TYPE> & ElemStoreSol<TYPE>::operator= (const ElemStoreSol & x)
{
  ENTER_FCN("ElemStoreSol::operator=(const ElemStoreSol &");
  Error("Not implemented here", __FILE__,__LINE__); 
}

template<class TYPE>
BaseStoreSol & ElemStoreSol<TYPE>::operator= (const BaseStoreSol & x)
{
  ENTER_FCN("ElemStoreSol::operator=(const ElemStoreSol &");
  if ( &x == dynamic_cast<ElemStoreSol*>(this))
    return (*this);

  const ElemStoreSol<TYPE> & temp = dynamic_cast<const ElemStoreSol<TYPE>&>(x);
  
  this->numNodes_ = temp.numNodes_;
  this->numSolutions_ = temp.numSolutions_;
  this->length_ = temp.length_;
  this->solTypes_ = temp.solTypes_;
  this->solOffset_ = temp.solOffset_;
  this->solDofs_ = temp.solDofs_;
  this->totalDofs_ = temp.totalDofs_;
  this->data_ = temp.data_;
  this->convertedData_ = temp.convertedData_;
 
  return dynamic_cast<ElemStoreSol &> (*this);
}

template <class TYPE>
void ElemStoreSol<TYPE>::Print(std::ostream& str)			  
{
  ENTER_IFCN( "operator<<(ElemStoreSol<TYPE>) ");
  Error( "not implemented" );

  std::map<SolutionType,Integer>::const_iterator it;
  //for (it = solDofs_.begin(); it!=solDofs_.end(); it++)
    
}





} //namespace
