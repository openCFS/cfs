#include "nodestoresol.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"

#include <typeinfo>

namespace CoupledField{



template<class TYPE>
NodeStoreSol<TYPE>::NodeStoreSol()
{
  ENTER_FCN( "NodeStoreSol::NOdeStoreSol()" );
  
  numNodes_ = 0;
  numSolutions_ = 0;
  length_ = 0;
  totalDofs_ = 0;
  eqnOffset_ = 0;

}
  
template<class TYPE>
NodeStoreSol<TYPE>::NodeStoreSol(Integer numNodes, 
			 StdVector<SolutionType> solTypes, 
			 StdVector<Integer> solDofs)
{
  ENTER_FCN( "NodeStoreSol::NodeStoreSol(numNodes, solTypes, solDofs" );
  Error( "Not implemented here", __FILE__, __LINE__ );
}

template<class TYPE>
NodeStoreSol<TYPE>::NodeStoreSol(const Integer numNodes,
			 const SolutionType solType,
			 const Integer numDofs)
{
  ENTER_FCN( "NodeStoreSol::NodeStoreSol(numNodes, solType, soDofs" );
  Error( "Not implemented here", __FILE__, __LINE__ );
}


template<class TYPE>
NodeStoreSol<TYPE>::NodeStoreSol(const NodeStoreSol & x) 
{
  ENTER_FCN( "NodeStoreSol::NodeStoreSol(const NodeStoreSol)" );
  Error( "Not implemented here", __FILE__, __LINE__ );
}



template<class TYPE>
NodeStoreSol<TYPE>::~NodeStoreSol() 
{
  ENTER_FCN( "NodeStoreSol::~NodeStoreSol" );
 
  
}

template<class TYPE>
void NodeStoreSol<TYPE>::Clear()
{
  ENTER_FCN( "NodeStoreSol::Clear()" );
}

template<class TYPE>
void NodeStoreSol<TYPE>::Init()  
{
  Init(TYPE());
}

template<class TYPE>
void NodeStoreSol<TYPE>::SetPtrEQNData(NodeEQN * ptNodeEQN)
{
  ENTER_FCN( "NodeStoreSol::SetPtrEQNData" );
  ptEQN_ = ptNodeEQN;
}

template<class TYPE>
void NodeStoreSol<TYPE>::Init(const TYPE val)  
{
  ENTER_FCN( "NodeStoreSol::Init" );

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
      Error("NodeStoreSol::Init(): Before calling Init(), the number of solutions,\
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

      // Determine the 'real' length of
      // the solution vector, which depends on the 
      // type of equation mapping
      if (ptEQN_->IsBlockMapped())
	{
	  lengthVector_ = ptEQN_->GetNumEQNs() * totalDofs_;
	  eqnDofs_ = totalDofs_;
	} 
      else {
	lengthVector_ = ptEQN_->GetNumEQNs();
	eqnDofs_ = 1;
      }
   
      data_.Resize(lengthVector_);
    }

  data_.Init(val);
}

template<class TYPE>
void NodeStoreSol<TYPE>::SetNumSolutions(const Integer nSols)
{
  ENTER_IFCN("NodeStoreSol::SetNumSolutions");

  numSolutions_ = nSols;
  length_ = 0;
}

template<class TYPE>
void NodeStoreSol<TYPE>::SetNumNodes(const Integer nNodes)
{
  ENTER_IFCN("NodeStoreSol::SetNumNodes");
  numNodes_ = nNodes;
  length_ = 0;
}

template<class TYPE>
void NodeStoreSol<TYPE>::SetSolutionType(const SolutionType solType, const Integer numSolution)
{
  ENTER_IFCN("NodeStoreSol::SetSolutionType");

#ifdef CHECK_INDEX
  if (numSolution >= numSolutions_) Error("NodeStoreSol: Index out of Bounds",__FILE__,__LINE__);
#endif

  solTypes_[solType] = numSolution;
  length_ = 0;
}

template<class TYPE>
void NodeStoreSol<TYPE>::SetNumDofs(const Integer dof, const SolutionType sol)
{
  ENTER_IFCN("NodeStoreSol::SetDof");
  
  // check if only one dof was assigned
  // -> only one entry exists
  if (numSolutions_ == 1 && sol == NO_SOLUTION_TYPE)
    solDofs_[(*solTypes_.begin()).first] = dof;
  else
    solDofs_[sol] = dof;

  length_ = 0;
}

template<class TYPE>
void NodeStoreSol<TYPE>::GetSolutionTypes(StdVector<SolutionType> &solTypes) const
{
  ENTER_FCN( "NodeStoreSol::GetSolutionTypes" );
  
  solTypes.Resize(numSolutions_);
  std::map<SolutionType,Integer>::const_iterator it;
 
    for (it = solTypes_.begin(); it!=solTypes_.end(); it++)
      solTypes[(*it).second] = (*it).first;
  
}

template<class TYPE>
TYPE& NodeStoreSol<TYPE>::operator()(Integer node, Integer dof)
{
  ENTER_IFCN("NodeStoreSol::operator()");
 #ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (numSolutions_ > 1) 
    Error("NodeStoreSol:operator(): Only defined objects with one type of solution!",__FILE__,__LINE__);
#endif
  Integer eqnNr,eqnDof;  
  ptEQN_->Node2EQN(node+1,dof+1,eqnNr,eqnDof);

  Warning ("Is this operator ever be used?", __FILE__, __LINE__);
  if (eqnNr > 0)
    return data_[abs(eqnNr-1)+eqnDof-1];
  else
    Error("NodeStoreSol::operator(): This operator gives only writing access to non-BC nodes.",
	  __FILE__, __LINE__);
}

template<class TYPE>
TYPE  NodeStoreSol<TYPE>:: operator()(Integer node, Integer dof) const
{
  ENTER_IFCN("NodeStoreSol::operator()");
  
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (numSolutions_ > 1) 
    Error("NodeStoreSol:operator(): Only defined objects with one type of solution!",__FILE__,__LINE__);
#endif
  
  Integer eqnNr,eqnDof;  
  ptEQN_->Node2EQN(node+1,dof+1,eqnNr,eqnDof);

  Warning ("Is this operator ever be used?", __FILE__, __LINE__);
  if (eqnNr > 0)
    return data_[abs(eqnNr-1)+eqnDof-1];
  else
    return TYPE();
}


template<class TYPE>
void NodeStoreSol<TYPE>::GetSolVector(const SolutionType type, CFSVector & val) const
{
  ENTER_FCN("NodeStoreSol::GetSolVector");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
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
void NodeStoreSol<TYPE>::SetSolVector(const SolutionType type, const CFSVector & val)
{
  ENTER_FCN("NodeStoreSol::SetSolVector");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif 
  
#ifdef CHECK_INDEX
  if (val.GetSize() != data_.GetSize())
    Error("NodeStoreSol::SetSolVector(): Vector has incompatible dimensions!",__FILE__,__LINE__);
#endif

  const Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
  for (Integer i=0; i<temp.GetSize(); i++)
    data_[i] = temp[i];
}

template<class TYPE>
void NodeStoreSol<TYPE>::GetSolution(const SolutionType type, BaseStoreSol & val) const
{
  ENTER_FCN("NodeStoreSol::GetSolution");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Integer offset = (*solOffset_.find(type)).second;
  Integer dof = (*solDofs_.find(type)).second;
    
  NodeStoreSol<TYPE> & temp = dynamic_cast<NodeStoreSol<TYPE>&>(val);

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
  temp.eqnOffset_ = offset;
  temp.ptEQN_ = ptEQN_;
  temp.length_ = temp.numNodes_ * dof;
  temp.data_.Resize(temp.length_);

  // TEMPORARYLY
  Integer myLength = (Integer) lengthVector_ / totalDofs_;


  for (Integer iNode=0; iNode<myLength; iNode++)
    for (Integer iDof=0; iDof<dof; iDof++)
      temp.data_[iNode*dof+iDof] = data_[iNode*totalDofs_+iDof+offset];
}

template<class TYPE>
void NodeStoreSol<TYPE>::SetNodalResult(const Integer nodeNr, const CFSVector &val)
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
  for (Integer i=0; i<temp.GetSize(); i++)
    data_[nodeNr*totalDofs_ + i] = temp[i];
}

template<class TYPE>
void NodeStoreSol<TYPE>::GetNodalResult(const Integer nodeNr, CFSVector &val) const
{
  ENTER_FCN("NodeStoreSol::GetNodalResult");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Error("Not implemented here", __FILE__,__LINE__);
}

template<class TYPE>
void NodeStoreSol<TYPE>::GetSolVectorSingleDof(const SolutionType type, const Integer dof, CFSVector & val) const
{
  ENTER_FCN("NodeStoreSol::GetSolVectorSingleDof");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Error( "Not implemented here" );
}

template<class TYPE>
void NodeStoreSol<TYPE>::GetSolVectorSingleDof(const Integer dof, CFSVector & val) const
{
  ENTER_FCN("NodeStoreSol::GetSolVectorSingleDof");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
  Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
 
  Integer globalPos, eqnDofs;
  
  temp.Resize(ptEQN_->GetNumGlobalNodes());
  eqnDofs = ptEQN_->GetNumDofs();
 
  Warning ("Not sure if this is working properly!!!");
  // Loop over all Equations
  for (Integer iEQN=0; iEQN<ptEQN_->GetNumEQNs(); iEQN++)
    {
      ptEQN_->EQN2SolVectorPos(iEQN+1, 1, globalPos);

      // ONLY TEMPORARY, until superblocknumbering
      // for piezoPDE works ;-)
      globalPos =(Integer) globalPos/eqnDofs;
      
      temp[globalPos] = data_[iEQN*totalDofs_ + dof];
    }
}


template<class TYPE>
void NodeStoreSol<TYPE>::Get(const Integer nodeNr, const Integer dof, TYPE & ret) const
{
  ENTER_FCN("NodeStoreSol::Get");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (numSolutions_ > 1)
    Error("NodeStoreSol::Get(): Only used for single solution objects!",__FILE__,__LINE__);
  if (nodeNr > numNodes_)
    Error("NodeStoreSol::Get(): Index out of bounds",__FILE__,__LINE__);
#endif

  Integer eqnNr, eqnDof;
  ptEQN_->Node2EQN(nodeNr+1,dof+1,eqnNr,eqnDof);

  if (eqnNr != 0)
    ret = data_[abs(eqnNr-1)+eqnDof-1];
  else
    ret =  TYPE();
}

template<class TYPE>
void NodeStoreSol<TYPE>::Get(const SolutionType type, const Integer nodeNr, const Integer dof, TYPE & ret) const
{
  ENTER_FCN("NodeStoreSol::Get");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Error ("Not yet adapted to EQN class", __FILE__, __LINE__);
  Integer offset = (*solOffset_.find(type)).second;

  ret = data_[nodeNr * totalDofs_ + offset + dof];
}

template<class TYPE>
void NodeStoreSol<TYPE>::Set(const SolutionType type, const Integer nodeNr, const Integer dof, const TYPE val)
{
  ENTER_FCN("NodeStoreSol::Set");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Error("Not implemented here", __FILE__,__LINE__);
}

template<class TYPE>
void NodeStoreSol<TYPE>::Add(const SolutionType type, const Integer nodeNr, const Integer dof, const TYPE val) const
{
  ENTER_FCN("NodeStoreSol::Set");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  Error("Not implemented here", __FILE__,__LINE__);
}

template<class TYPE>
void NodeStoreSol<TYPE>::SetCompleteVector(const CFSVector & val)
{
  ENTER_FCN("NodeStoreSol::SetCompleteVector");
#ifdef CHECK_INITIALIZED
   if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
#ifdef CHECK_INDEX
   if (val.GetSize() !=  lengthVector_)
     Error("NodeStoreSol::SetCompleteVector(): Vector has wrong size!",__FILE__,__LINE__);
#endif

   const  Vector<TYPE> & temp = dynamic_cast<const Vector<TYPE>&>(val);
   data_ = temp;
}
  
template<class TYPE> 
void NodeStoreSol<TYPE>::GetCompleteVector(CFSVector & val) const
{
 ENTER_FCN("NodeStoreSol::SetCompleteVector");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
  Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(val);
  temp = data_;
}

template<class TYPE>
void NodeStoreSol<TYPE>::GetVectorPointer(CFSVector* &ptrToVec)
{
  ENTER_FCN("NodeStoreSol::GetVectorPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  ptrToVec = (CFSVector*) &data_;
}

template<class TYPE>
void NodeStoreSol<TYPE>::GetVectorPointer(Vector<TYPE>* &ptrToVec)
{
  ENTER_FCN("NodeStoreSol::GetVectorPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  ptrToVec =  &data_;
}

template<class TYPE>
void NodeStoreSol<TYPE>::CopyFromDataPointer(TYPE * ptr)
{
  ENTER_FCN("NodeStoreSol::CopyFromDataPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  for (Integer i=0; i<lengthVector_; i++)
    {
      data_[i] = ptr[i];
    }

}

template<class TYPE>
void NodeStoreSol<TYPE>::SetDataPointer(TYPE * ptr)
{
  ENTER_FCN("NodeStoreSol::SetDataPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",
				__FILE__,__LINE__);
#endif
  
  data_.data_ = ptr;

}

template<class TYPE>
void NodeStoreSol<TYPE>::GetDataPointer(TYPE* &ptr)
{
  ENTER_FCN("NodeStoreSol::GetDataPointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",
				__FILE__,__LINE__);
#endif
 
  ptr =  data_.data_;
}


template<>
Double* NodeStoreSol<Double>::GetDoublePointer()
{
  ENTER_FCN("NodeStoreSol::GetDoublePointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",
				__FILE__,__LINE__);
#endif
  return data_.data_;
}


template<class TYPE>
Double* NodeStoreSol<TYPE>::GetDoublePointer()
{
  ENTER_FCN("NodeStoreSol::GetDoublePointer");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  Error("Not implemented here",__FILE__,__LINE__);
}



///////// Transformation Operations ///////// 

template<class TYPE>
void NodeStoreSol<TYPE>::GetElemSolution(CFSVector & elemSol,
			       const StdVector<Integer> & connect) const
{
  ENTER_FCN( "NodeStoreSol::GetElemSolution" );
  Integer eqnNr, eqnDof;

  Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(elemSol);
  
  temp.Resize(totalDofs_*connect.GetSize());

  for (Integer iDof=0; iDof<totalDofs_; iDof++)
    for (Integer iNode=0; iNode<connect.GetSize(); iNode++)
      {
	ptEQN_->Node2EQN(connect[iNode],iDof+1,eqnNr,eqnDof);
	if (eqnNr != 0)
	  //temp[iDof][iNode] = data_[totalDofs_*(connect[iNode]-1) + iDof];
	  temp[iDof+iNode*totalDofs_] = data_[abs(eqnNr-1)*eqnDofs_ + eqnDof-1];
	else
	  temp[iDof + iNode*totalDofs_] = TYPE();
      }
}


template<class TYPE>
void NodeStoreSol<TYPE>::GetElemSolutionAsMatrix(CFSMatrix & elemSol, 
					     const StdVector<Integer> & connect) const
{
  ENTER_FCN("NodeStoreSol::GetElemSolutionAsMatrix");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif
  
  Integer eqnNr, eqnDof;

  Matrix<TYPE> & temp = dynamic_cast<Matrix<TYPE>&>(elemSol);
  
  temp.Resize(totalDofs_,connect.GetSize());

  for (Integer iDof=0; iDof<totalDofs_; iDof++)
    for (Integer iNode=0; iNode<connect.GetSize(); iNode++)
      {
	ptEQN_->Node2EQN(connect[iNode],iDof+1,eqnNr,eqnDof);
	if (eqnNr != 0)
	  temp[iDof][iNode] = data_[abs(eqnNr-1)*eqnDofs_ + eqnDof-1];
	else
	  temp[iDof][iNode] = TYPE();
      }
}

template<class TYPE>
void NodeStoreSol<TYPE>::TransformNodeSolution(CFSVector & transformedSolution,
					       Grid * ptGrid,
					       const Integer level) const
{
  ENTER_FCN("NodeStoreSol::TransformNodeSolution");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) 
    Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  //std::cerr << std::endl << std::endl;
  //  std::cerr << "Entering TransformNodeSolution" << std::endl;
  
  Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(transformedSolution);
 
  Integer globalPos, eqnDofs, factor, dofsPerEQN, eqnDof;
  
  temp.Resize(ptGrid->GetMaxnumnodes(level)*totalDofs_);
  eqnDofs = ptEQN_->GetNumDofs();
  dofsPerEQN = ptEQN_->GetNumDofsPerEQN();
 
  //std::cerr << "temp has size " << temp.GetSize() << std::endl;
  //std::cerr << "data_ has size " << data_.GetSize() << std::endl;
  //std::cerr << "totalDofs = " << totalDofs_ << std::endl;


   if (ptEQN_->IsBlockMapped())
     {
       // In this case each eqn has a dof number > 1 and we need
       // two loops to map each eqn to its positionS
       
       // Loop over all Equations
       for (Integer iEQN=0; iEQN<(Integer) ptEQN_->GetNumEQNs(); iEQN++)
	 for (Integer iDof=0; iDof<totalDofs_; iDof++)
	   {
	     ptEQN_->EQN2SolVectorPos(iEQN+1,eqnOffset_+iDof+1,globalPos);
	     //std::cerr << "SolVectorPos for EQN " << iEQN+1 << " is " << globalPos << std::endl;
	     // ONLY TEMPORARY, until superblocknumbering
	     // for piezoPDE works ;-)
	     globalPos =(Integer) (globalPos-iDof)/eqnDofs *totalDofs_+iDof;    bool isDamping_ = false;

	  
	     //std::cerr << "Corrected position is " << globalPos << std::endl;
	     // 	  for (Integer iDof=0; iDof<totalDofs_; iDof++)
	     // 	    {
	     //std::cerr << "Equation" << iEQN+1 << "gets stored to location" << globalPos+iDof << std::endl;
	     //std::cerr << "Data is read from data_[" << iEQN*totalDofs_+iDof << "] = " <<  data_[iEQN*totalDofs_+iDof] << std::endl;
	     //std::cerr << "Try to store data_[" <<  +iDof;
	     
	     //    std::cerr <<  "] to temp[" << globalPos+iDof << "]" << std::endl;
	     temp[globalPos] = data_[iEQN*totalDofs_ +iDof];
	     // }
	   }
     } 
   else 
     {
       // In this case each eqn has only one dof 
       // and only one for loop is needed
      for (Integer iEQN=0; iEQN<(Integer) ptEQN_->GetNumEQNs(); iEQN++)
	  {
	  ptEQN_->EQN2SolVectorPos(iEQN+1,1,globalPos);
	  //std::cerr << "SolVectorPos for EQN " << iEQN+1 << " is " << globalPos << std::endl;
	  // ONLY TEMPORARY, until superblocknumbering
	  // for piezoPDE works ;-)
	  //globalPos =(Integer) globalPos/eqnDofs *totalDofs_;
	  
	  //std::cerr << "Corrected position is " << globalPos << std::endl;
	
	  temp[globalPos] = data_[iEQN];
	  // }
	}
      
    }
}  

template<class TYPE>
void NodeStoreSol<TYPE>::NodeSolutionToCoupling(CFSVector & couplingSol,
						const StdVector<Integer>& nodeNumbers) const
{
  ENTER_FCN("NodeStoreSol:::NodeSolutionToCoupling");
#ifdef CHECK_INITIALIZED
  if (length_ == 0) Error("NodeStoreSol: Use of uninitialized object!",__FILE__,__LINE__);
#endif

  Vector<TYPE> & temp = dynamic_cast<Vector<TYPE>&>(couplingSol);
  Integer eqnNr, eqnDof;

  temp.Resize(nodeNumbers.GetSize() * totalDofs_);
  for (Integer iNode=0; iNode<nodeNumbers.GetSize(); iNode++)
    for (Integer iDof=0; iDof<totalDofs_; iDof++)
      {
	ptEQN_->Node2EQN(nodeNumbers[iNode], iDof+1, eqnNr, eqnDof);
	if (eqnNr != 0)
	  temp.data_[iNode*totalDofs_ + iDof] = data_[abs((eqnNr-1)*eqnDofs_ + eqnDof-1)];
	else
	  temp.data_[iNode*totalDofs_ + iDof] = TYPE();
	
	//std::cerr << "TocOupling[" << nodeNumbers[iNode] << " = " << data_[abs((eqn-1)*totalDofs_+iDof)] << std::endl;
      }
}

  
//////////////////// Operators //////////////////
template<class TYPE>
NodeStoreSol<TYPE> & NodeStoreSol<TYPE>::operator= (const NodeStoreSol & x)
{
  ENTER_FCN("NodeStoreSol::operator=(const NodeStoreSol &");
  Error("Not implemented here", __FILE__,__LINE__); 
}

template<class TYPE>
BaseStoreSol & NodeStoreSol<TYPE>::operator= (const BaseStoreSol & x)
{
  ENTER_FCN("NodeStoreSol::operator=(const BaseNodeStoreSol &");
  if ( &x == dynamic_cast<BaseStoreSol*>(this))
    return (*this);

  const NodeStoreSol<TYPE> & temp = dynamic_cast<const NodeStoreSol<TYPE>&>(x);
  
  this->numNodes_ = temp.numNodes_;
  this->numSolutions_ = temp.numSolutions_;
  this->length_ = temp.length_;
  this->solTypes_ = temp.solTypes_;
  this->solOffset_ = temp.solOffset_;
  this->solDofs_ = temp.solDofs_;
  this->totalDofs_ = temp.totalDofs_;
  this->data_ = temp.data_;
  this->convertedData_ = temp.convertedData_;
 
  return dynamic_cast<BaseStoreSol &> (*this);
}

template <class TYPE>
void NodeStoreSol<TYPE>::Print(std::ostream& str)			  
{
  ENTER_IFCN( "operator<<(NodeStoreSol<TYPE>) ");
  Error( "not implemented" );

  std::map<SolutionType,Integer>::const_iterator it;
  
}





} //namespace
