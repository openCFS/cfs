#include "storesol.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField{



template<class TYPE>
StoreSol<TYPE>::StoreSol()
{
  ENTER_FCN("StoreSol::StoreSol()");
  
  numNodes_ = 0;
  numSolutions_ = 0;
  length_ = 0;
  totalDofs_ = 0;

}
  
template<class TYPE>
StoreSol<TYPE>::StoreSol(Integer numNodes, 
			 std::vector<SolutionType> solutionTypes, 
			 std::vector<Integer> solutionDofs)
{
  ENTER_FCN("StoreSol::StoreSol(numNodes, solutiontypes,solutionDofs");
  Info->Error("Not implemented here", __FILE__,__LINE__);
}

template<class TYPE>
StoreSol<TYPE>::StoreSol(const StoreSol & x) 
{
  ENTER_FCN("StoreSol::StoreSol(const StoreSol)");
  Info->Error("Not implemented here", __FILE__,__LINE__);
}



template<class TYPE>
StoreSol<TYPE>::~StoreSol() 
{
  ENTER_FCN("StoreSol::~StoreSol");
 
  
}


template<class TYPE>
void StoreSol<TYPE>::Init(const TYPE val)  
{
  ENTER_FCN("StoreSol::Init");
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
      Info->Error("StoreSol::Init(): Before calling Init(), the number of solutions,\
           nodes, types and dofs has to be set to NONZERO value!",__FILE__,__LINE__);
    }
  
#endif
  
  // only the first time the struct gets initialized
  if (length_ == 0)
    {
      if (solDofs_.size() != numSolutions_ || solTypes_.size() != numSolutions_)
	Info->Error("Inconsistent definition of Storesolution class",__FILE__,__LINE__);
      
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
void StoreSol<TYPE>::SetNumSolutions(const Integer nSols)
{
  ENTER_IFCN("StoreSol::SetNumSolutions");

  numSolutions_ = nSols;

  length_ = 0;
  solTypes_.clear();
  solDofs_.clear();
  solOffset_.clear();
}

template<class TYPE>
void StoreSol<TYPE>::SetNumNodes(const Integer nNodes)
{
  ENTER_IFCN("StoreSol::SetNumNodes");
  numNodes_ = nNodes;
}

template<class TYPE>
void StoreSol<TYPE>::SetSolutionType(const SolutionType solType, const Integer numSolution)
{
  ENTER_IFCN("StoreSol::SetSolutionType");

#ifdef CHECK_INDEX
  if (numSolution >= numSolutions_) Info->Error("StoreSol: Index out of Bounds",__FILE__,__LINE__);
#endif

  solTypes_[solType] = numSolution;
}

template<class TYPE>
void StoreSol<TYPE>::SetNumDofs(const Integer dof, const SolutionType sol)
{
  ENTER_IFCN("StoreSol::SetDof");
  
  // check if only one dof was assigned
  // -> only one entry exists
  if (numSolutions_ == 1 && sol == NO_SOLUTION_TYPE)
    solDofs_[(*solTypes_.begin()).first] = dof;
  else
    solDofs_[sol] = dof;
}

template<class TYPE>
TYPE& StoreSol<TYPE>:: operator()(Integer node, Integer dof)
{
  ENTER_IFCN("StoreSol::operator()");
  
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (numSolutions_ > 1) 
    Info->Error("StoreSol:operator(): Only defined objects with one type of solution!",__FILE__,__LINE__);
#endif
  return data_[node*totalDofs_ + dof];
}

template<class TYPE>
TYPE  StoreSol<TYPE>:: operator()(Integer node, Integer dof) const
{
  ENTER_IFCN("StoreSol::operator()");
  
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (numSolutions_ > 1) 
    Info->Error("StoreSol:operator(): Only defined objects with one type of solution!",__FILE__,__LINE__);
#endif
  return data_[node*totalDofs_ + dof];
}


template<class TYPE>
void StoreSol<TYPE>::GetSolVector(const SolutionType type, CFSVector & val) const
{
  ENTER_FCN("StoreSol::GetSolVector");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
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
void StoreSol<TYPE>::SetSolVector(const SolutionType type, const CFSVector & val)
{
  ENTER_FCN("StoreSol::SetSolVector");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif 
  
#ifdef CHECK_INDEX
  if (val.GetSize() != data_.GetSize())
    Info->Error("StoreSik:SetSolVector(): Incompatible dimensions",__FILE__,__LINE__);
#endif

  const Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
  for (Integer i=0; i<temp.GetSize(); i++)
    data_[i] = temp[i];
}

template<class TYPE>
void StoreSol<TYPE>::GetSolution(const SolutionType type, BaseStoreSol & val) const
{
  ENTER_FCN("StoreSol::GetSolution");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Integer offset = (*solOffset_.find(type)).second;
  Integer dof = (*solDofs_.find(type)).second;
    
  StoreSol<TYPE> & temp = dynamic_cast<StoreSol<TYPE>&>(val);

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
void StoreSol<TYPE>::SetNodalResult(const Integer nodeNr, const CFSVector &val)
{
  ENTER_FCN("StoreSol::SetNodalResult");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
#ifdef CHECK_INDEX
  if (nodeNr >= numNodes_) 
    Info->Error("StoreSol::SetNodalResult(): index out of bounds",__FILE__,__LINE__);
  if (val.GetSize() != totalDofs_)
    Info->Error("StoreSol::SetNodalResult(): vector of incompatible dimension",__FILE__,__LINE__);
#endif

  const Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
  for (Integer i=0; i<temp.GetSize(); i++)
    data_[nodeNr*totalDofs_ + i] = temp[i];
}

template<class TYPE>
void StoreSol<TYPE>::GetNodalResult(const Integer nodeNr, CFSVector &val) const
{
  ENTER_FCN("StoreSol::GetNodalResult");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Info->Error("Not implemented here", __FILE__,__LINE__);
}

template<class TYPE>
void StoreSol<TYPE>::GetSolVectorSingleDof(const SolutionType type, const Integer dof, CFSVector & val) const
{
  ENTER_FCN("StoreSol::GetSolVectorSingleDof");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Info->Error("Not implemented here", __FILE__,__LINE__);
}

template<class TYPE>
void StoreSol<TYPE>::GetSolVectorSingleDof(const Integer dof, CFSVector & val) const
{
  ENTER_FCN("StoreSol::GetSolVectorSingleDof");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Info->Error("Not implemented here", __FILE__,__LINE__);
}


template<class TYPE>
void StoreSol<TYPE>::Get(const Integer nodeNr, const Integer dof, TYPE & ret) const
{
  ENTER_FCN("StoreSol::Get");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (numSolutions_ > 1)
    Info->Error("StoreSol::Get(): Only used for single solution objects!",__FILE__,__LINE__);
  if (nodeNr > numNodes_)
    Info->Error("StoreSol::Get(): Index out of bounds",__FILE__,__LINE__);
#endif

  ret = data_[totalDofs_*nodeNr + dof];
}

template<class TYPE>
void StoreSol<TYPE>::Get(const SolutionType type, const Integer nodeNr, const Integer dof, TYPE & ret) const
{
  ENTER_FCN("StoreSol::Get");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Integer offset = (*solOffset_.find(type)).second;

  ret = data_[nodeNr * totalDofs_ + offset + dof];
}

template<class TYPE>
void StoreSol<TYPE>::Set(const SolutionType type, const Integer nodeNr, const Integer dof, const TYPE val)
{
  ENTER_FCN("StoreSol::Set");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Info->Error("Not implemented here", __FILE__,__LINE__);
}

template<class TYPE>
void StoreSol<TYPE>::Add(const SolutionType type, const Integer nodeNr, const Integer dof, const TYPE val) const
{
  ENTER_FCN("StoreSol::Set");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Info->Error("Not implemented here", __FILE__,__LINE__);
}

template<class TYPE>
void StoreSol<TYPE>::SetCompleteVector(const CFSVector & val)
{
  ENTER_FCN("StoreSol::SetCompleteVector");
#ifdef CHECK_INITIALIZED
   if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
#ifdef CHECK_INDEX
   if (val.GetSize() !=  length_)
     Info->Error("StoreSol::SetCompleteVector(): Vector has wrong size!",__FILE__,__LINE__);
#endif

   const  Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
   data_ = temp;
}
  
template<class TYPE> 
void StoreSol<TYPE>::GetCompleteVector(CFSVector & val) const
{
 ENTER_FCN("StoreSol::SetCompleteVector");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
  Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
  temp = data_;
}

template<class TYPE>
void StoreSol<TYPE>::GetVectorPointer(CFSVector* &ptrToVec)
{
  ENTER_FCN("StoreSol::GetVectorPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  ptrToVec = (CFSVector*) &data_;
}

template<class TYPE>
void StoreSol<TYPE>::GetVectorPointer(Vector<TYPE>* &ptrToVec)
{
  ENTER_FCN("StoreSol::GetVectorPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  ptrToVec =  &data_;
}

template<class TYPE>
void StoreSol<TYPE>::CopyFromDataPointer(TYPE * ptr)
{
  ENTER_FCN("StoreSol::CopyFromDataPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  for (Integer i=0; i<length_; i++)
    data_[i] = ptr[i];

}

template<class TYPE>
void StoreSol<TYPE>::SetDataPointer(TYPE * ptr)
{
  ENTER_FCN("StoreSol::SetDataPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
  // Unsicher ...
  //if(data_.data_)
  // delete[] data_.data_;
 
  data_.data_ = ptr;

}

template<class TYPE>
void StoreSol<TYPE>::GetDataPointer(TYPE* &ptr)
{
  ENTER_FCN("StoreSol::GetDataPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
 
  ptr =  data_.data_;
}


template<>
Double* StoreSol<Double>::GetDoublePointer()
{
  ENTER_FCN("StoreSol::GetDoublePointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  return data_.data_;
}


template<class TYPE>
Double* StoreSol<TYPE>::GetDoublePointer()
{
  ENTER_FCN("StoreSol::GetDoublePointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  Info->Error("Not implemented here",__FILE__,__LINE__);
}



///////// Transformation Operations ///////// 
template<class TYPE>
void StoreSol<TYPE>::GetElemSolutionAsMatrix(CFSMatrix & elemSol, Vector<Integer> & connect) const
{
  ENTER_FCN("StoreSol::GetElemSolutionAsMatrix");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
  Matrix<TYPE> & temp = dynamic_cast<Matrix<TYPE>&>(elemSol);
  temp.Resize(totalDofs_,connect.GetSize());
  for (Integer iDof=0; iDof<totalDofs_; iDof++)
    for (Integer iNode=0; iNode<connect.GetSize(); iNode++)
      temp[iDof][iNode] = data_[totalDofs_*(connect[iNode]-1) + iDof];
}

template<class TYPE>
void StoreSol<TYPE>::TransformNodeSolution(BaseStoreSol & transformedSolution,
					   const std::vector<Integer> & mapping,
					   Grid * ptGrid,
					   const Integer level) const
{
  ENTER_FCN("StoreSol::TransformNodeSolution");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) 
    Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  StoreSol<TYPE> & temp = dynamic_cast<StoreSol<TYPE>&>(transformedSolution);
 
  // initialize the return basesolution
  //std::cerr << "numNodes: " << numNodes_ << std::endl;
  //std::cerr << "solDofs: " << solDofs_ << std::endl;
  //std::cerr << "numSolutions: " << numSolutions_ << std::endl;
  //std::cerr << "totalDofs: " << totalDofs_ << std::endl;
  //std::cerr << "length: " << length_ << std::endl;
  

  temp.numNodes_ = ptGrid->GetMaxnumnodes(level);
  temp.solDofs_ = solDofs_;
  temp.solTypes_ = solTypes_;
  temp.solOffset_ = solOffset_;
  temp.numSolutions_ = numSolutions_;
  temp.totalDofs_ = totalDofs_;
  temp.length_ = temp.numNodes_ * totalDofs_;
  temp.data_.Resize(temp.length_);

  // Loop over all PDE nodes
  for (Integer iNode=0; iNode<mapping.size(); iNode++)
    // loop over dimensions
     for (Integer iDof=0; iDof<totalDofs_; iDof++)
       temp.data_[(mapping[iNode]-1)*totalDofs_ + iDof] = data_[iNode*totalDofs_ + iDof];
}  

template<class TYPE>
void StoreSol<TYPE>::TransformElemSolution(BaseStoreSol & transformedSolution,
					   const std::vector<Integer> & mapping,
					   Grid * ptGrid,
					   const Integer level) const
{
  ENTER_FCN("StoreSol::TransformElemSolution");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
  StoreSol<TYPE> & temp = dynamic_cast<StoreSol<TYPE>&>(transformedSolution);
  temp.numNodes_ = ptGrid->GetMaxnumElem(level);
  temp.solDofs_ = solDofs_;
  temp.solTypes_ = solTypes_;
  temp.solOffset_ = solOffset_;
  temp.numSolutions_ = numSolutions_;
  temp.totalDofs_ = totalDofs_;
  temp.length_ = temp.numNodes_ * totalDofs_;
  temp.data_.Resize(temp.length_); 

  // Loop over all PDE elements
  for (Integer iElem=0; iElem<mapping.size(); iElem++)
    // Loop over all dimensions
    for (Integer iDof=0; iDof<totalDofs_; iDof++)
      temp.data_[(mapping[iElem]-1)*totalDofs_ + iDof] = data_[iElem*totalDofs_ + iDof];
}

template<class TYPE>
void StoreSol<TYPE>::NodeSolutionToCoupling(BaseStoreSol & couplingSol,
					    const std::vector<Integer>& nodeNumbers,
					    const std::vector<Integer> & mapping) const
{
  ENTER_FCN("StoreSol:::NodeSolutionToCoupling");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  StoreSol<TYPE> & temp = dynamic_cast<StoreSol<TYPE>&>(couplingSol);
  
  temp.numNodes_ = nodeNumbers.size();
  temp.solDofs_ = solDofs_;
  temp.solTypes_ = solTypes_;
  temp.solOffset_ = solOffset_;
  temp.numSolutions_ = numSolutions_;
  temp.totalDofs_ = totalDofs_;
  temp.length_ = temp.numNodes_ * totalDofs_;
  temp.data_.Resize(temp.length_);

  for (Integer iDof=0; iDof<temp.totalDofs_; iDof++)
     for (Integer iNode=0; iNode<nodeNumbers.size(); iNode++)
       temp.data_[iNode*totalDofs_ + iDof] = data_[(mapping[nodeNumbers[iNode]-1 ] - 1)*totalDofs_ + iDof];
}

template<class TYPE> 
void StoreSol<TYPE>::ElemSolutionToCoupling(BaseStoreSol & couplingSol,
					    const std::vector<Elem*>& elements,
					    const CFSVector & elemSol) const
{
  ENTER_FCN("StoreSol::ElemSolutionToCoupling");

#ifdef CHECK_INITIALIZED
  if (length_ == 0) Info->Error("StoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  Info->Error("Not implemented here", __FILE__,__LINE__); 

}
  
  
//////////////////// Operators //////////////////
template<class TYPE>
StoreSol<TYPE> & StoreSol<TYPE>::operator= (const StoreSol & x)
{
  ENTER_FCN("StoreSol::operator=(const StoreSol &");
  Info->Error("Not implemented here", __FILE__,__LINE__); 
}

template<class TYPE>
BaseStoreSol & StoreSol<TYPE>::operator= (const BaseStoreSol & x)
{
  ENTER_FCN("StoreSol::operator=(const BaseStoreSol &");
}

} //namespace
