// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BASENODESTORESOL_2004
#define FILE_BASENODESTORESOL_2004

#include <iostream>
#include <map>

#include "Utils/StdVector.hh"
#include "DataInOut/WriteInfo.hh"
#include "General/environment.hh"
#include "Utils/tools.hh"
#include "Utils/vector.hh"
#include "PDE/eqnMap.hh"


namespace CoupledField{

  // Forward class declarations
  class Elem;
  class Grid;
  class CFSVector;
  class CFSMatrix;
  template<class TYPE> class Matrix;
  template<class TYPE> class Vector;


  //! This class is the new interface for handling solutions of PDEs

  //! This class is the new interface for handling the 
  //! (multidimensional) solutions of PDEs (instead of 
  //! old Array-class).
  //!
  //! In principle this class is only a wrapper around
  //! a CFSVector class. Additionally it contains
  //! informations about the information stored in it and
  //! it can selectively be accessed, set and read.
  //!
  //! A StoreSolution example for the Piezoelectric PDE would look
  //! like this:
  //! \f[ \left( \begin{array}{c} disp_x^1 \\ disp_y^1 \\  disp_z^1 \\  
  //! V_{elec}^1 \\ disp_x^2 \\ disp_y^2 \\  disp_z^2 \\ V_{elec}^2 \\ \cdots
  //! \end{array} \right) \f]
  //! Here the first three entries are the displacements for the first node,
  //! afterwards follows the electric Potential of it. This repeats for each
  //! node, so this is basically the layout of the solution vector as it is
  //! delivererd by the algebraic system. <br>
  //! The defining code for the above example would look like this:
  //! \verbatim
  //! mySol.SetNumSolutions(2); 
  //! mySol.SetSolutionType(MECH_DISPLACEMENT,0);
  //! mySol.SetSolutionType(ELEC_POTENTIAL,1);
  //! mySol.SetNumDofs(dof,MECH_DISPLACEMENT);
  //! mySol.SetNumDofs(1,ELEC_POTENTIAL);
  //! mySol.SetNumNodes(numNodes);
  //! mySol.Init(0.0) 
  //! \endverbatim
  //! \note An object of the StoreSolution can only used after the Init()-
  //! routine was called, otherwise an error is reported!
  //! \note Although the names of some methods refer to nodes, this class 
  //! also can handle element solutions.
  class BaseNodeStoreSol
  {
  public:
  
    //! Default destructor
  
    //! Since this is a base class used for inheritance we give it a virtual
    //! destructor.
    virtual ~BaseNodeStoreSol() {};
  
    //! Assignment operator

    //! Assignment operator has to be reimplemented in derived classes
    //! by the help of type casts.
    virtual BaseNodeStoreSol & operator= (const BaseNodeStoreSol & x) = 0;

    //! Set Pointer to nodal equation object
    virtual void SetPtrEQNData( EqnMap * eqnMap,
                                Grid *ptGrid) = 0;
    
    virtual void SetResult( shared_ptr<ResultInfo>& result ) {
      result_ = result; }
    
    virtual void SetRegions( StdVector<RegionIdType> regions ) {
      for( UInt i=0; i<regions.GetSize(); i++ ) {
        elems_.Push_back(shared_ptr<ElemList> (new ElemList( ptGrid_ ) ) );
        elems_.Last()->SetRegion( regions[i] );
      }
    }
          
    //! Hard coded query if values are complex
    virtual bool IsComplex() = 0;

    //! Deletes all data and layout information

    //! Deletes all data and layout information.
    //! This method should be called before the layout of the
    //! solution object is modified via SetNumSolution(), SetNumNodes(),
    //! SetNumDofs() or SetSolutionType().
    virtual void Clear() = 0;
  
    //! Initialization of the StoreSolution-object with 0-entries (REQUIRED)
  
    //! Initializes the object AFTER the other layout
    //! functions have been called with a given value
    //! \note Only after calling Init(), the object can
    //! store information
    virtual void Init() = 0;
  
    //! Initialization of the StoreSolution-object (REQUIRED)
  
    //! Initializes the object AFTER the other layout
    //! functions have been called with a given value
    //! \note Only after calling Init(), the object can
    //! store information
    /*!
      \param val (input) Value the object gets initialized with
    */
    virtual void Init(const Double val)
    {Error("BaseNodeStoreSol::Init() not implemented here", __FILE__, __LINE__);}
  
    //  virtual void Init(const Complex val)
    //  {Error("BaseNodeStoreSol::Init() not implemented here", __FILE__, __LINE__);}

    //! Set the number of different solution types
    /*! 
      \param nSols (input) Number of different solution types
    */
    //! \note All entries of this object are deleted afterwards
    virtual void SetNumSolutions(const UInt nSols) = 0;
  

    //! Set the number of solution nodes
    /*!
      \param nNodes (input) Number of solution nodes
    */
    //! \note All entries of this object are deleted afterwards
    virtual void SetNumNodes(const UInt nNodes) = 0;

  
    //! Set the type of solution and its position 

    //! Set the tpye of solution and its order of occurence
    //! in the solution vector (e.g. in Piezo, the mechanics displacement
    //! is stored previous to the electric potential)
    /*!
      \param solType (input) Type of solution (ref. enum SolutionType)
      \param numSolution (input) Occurence order in the solution vector
    */
    //! \note All entries of this object are deleted afterwards
    virtual void SetSolutionType(const SolutionType solType,
                                 UInt numSolution = 0) = 0;

  
    //! Set the number of dofs for one solution type

    //! Set the number of dofs for one solution.
    //! If only one type of solution is stored, the second parameter
    //! can be omitted.
    /*!
      \param dof (input) Number of degrees of freedom
      \param solType (input) Type of solution (ref. enum SolutionType)
    */
    //! \note All entries of this object are deleted afterwards
    virtual void SetNumDofs(const UInt dof,
                            const SolutionType solType = NO_SOLUTION_TYPE) = 0;


    //! Get list of stored solution types
    /*!
      \param solTypes (output) Vector of solution types
    */
    virtual void GetSolutionTypes(StdVector<SolutionType> &solTypes) const = 0;


    //! Return total number of dofs

    //! Return the accumulated number of Dofs for
    //! all solution types in this object
    inline UInt GetTotalNumDofs() const;

  
    //! Return the number of dofs for agiven solution type
  
    //! Return the number of dofs for a given solution type.
    //! If only one type of solution is stored, the second parameter
    //! can be omitted.
    /*!
      \param solType (input) Type of solution (ref. enum SolutionType)
    */
    inline UInt GetDof(const SolutionType solType = NO_SOLUTION_TYPE) const;

  
    //! Get number of nodes
    inline UInt GetNumNodes() const;
  

    //! Get total length of solution vector

    //! Get total length of solution vector (= totalDofs_ * numNodes_)
    inline UInt GetSize() const;

  
    //! Get number of solution types
    inline UInt GetNumSolutions() const;

  
    //! Get vectot with one type of solution. The length
    //! of the vector will be (numer of nodes in mesh *
    //! degree of freedoms)
    //! 
    /*!
      \param (input) Solution type (ref. enum SolutionType)
      \param (output) Vector with given solution type)
    */
    virtual void GetGlobalSolVector(const SolutionType solType, 
                                    CFSVector & val) const = 0;
  

    //! Set all solution types for one node
    /*!
      \param nodeNr (input) Node number for solution
      \param val (input) Vector containing nodal results
    */
    virtual void SetNodalResult(const UInt nodeNr,
                                const CFSVector &val) = 0;

  
    //! Get all solution types for one node
    /*!
      \param nodeNr (input) Node number for solution
      \param val (output) Vector containing nodal results
    */
    virtual void GetNodalResult(const UInt nodeNr,
                                CFSVector & val) const = 0;
  

    //! Get vector of one solution type for all nodes of one given dof
    /*!
      \param solType (input) Solution type (ref. enum SolutionType)
      \param dof (input)  Dof of solType
      \param val (output) Vector containing rsults
    */
    virtual void GetGlobalSolVectorSingleDof(const SolutionType solType,
                                             const UInt dof,
                                             CFSVector & val) const = 0;
                                         
  
  
    //! Get solution vector for all nodes  of one given dof
    /*!
      \param dof (input)  Dof of solType
      \param val (output) Vector containing rsults
    */
    //! \note This method may only be called if object contains only
    //! one type of solution.
    virtual void GetGlobalSolVectorSingleDof(const UInt dof,
                                             CFSVector & val) const = 0;


    //! Get single result of given node for given dof
    /*!
      \param nodeNr (input) Node number of result
      \param dof (input) Dof of result
      \param val (output) Result of node nodeNr for given dof
    */
    //! \note This method may only be called if object contains only
    //! one type of solution.
    virtual void Get(const UInt nodeNr, 
                     const UInt dof,
                     Double & val) const
    {Error("BaseNodeStoreSol::Get() not implemented here", __FILE__, __LINE__);} 

  
    //! Get single result of given solution type, node and dof
    /*! 
      \param solType (input) Solution type (ref. enum SolutionType)
      \param nodeNr (input) Node number of result
      \param dof (input) Dof of result
      \param val (output)  Result of node nodeNr for given dof
    */
    virtual void Get(const SolutionType solType, 
                     const UInt nodeNr,
                     const UInt dof, 
                     Double & val) const
    {Error("BaseNodeStoreSol::Get() not implemented here", __FILE__, __LINE__);}

  
    //! Set a single entry of a given solution type and a given dof
    /*!
      \param solType (input) Solution type (ref. enum SolutionType)
      \param nodeNr (input) Node number of result
      \param dof (input) Dof of result
      \param val (input)  Result of node nodeNr for given solType and dof
    */  
    virtual void Set(const SolutionType solType, 
                     const UInt nodeNr,
                     const UInt dof, 
                     const Double val)
    {Error("BaseNodeStoreSol::Set() not implemented here", __FILE__, __LINE__);}
  

    //! Add value to a single entry of a given solution type and a given dof
    /*!
      \param solType (input) Solution type (ref. enum SolutionType)
      \param nodeNr (input) Node number of result
      \param dof (input) Dof of result
      \param val (input)  Result of node nodeNr for given solType and dof
    */   
    virtual void Add(const SolutionType solType, 
                     const UInt nodeNr,
                     const UInt dof, 
                     const Double val) const
    {Error("BaseNodeStoreSol::Add() not implemented here", __FILE__, __LINE__);} 

  
    //! Set the complete solution vector inside this object
    /*!
      \param val (input) Vector containing complete solution
    */
    //! \note The layout of the val-vector must match the layout
    //! of the solution prescribed by SetNumSolutions(), SetNumDofs() and
    //! SetSolutionType()!
    virtual void SetAlgSysVector(const CFSVector & val) = 0;


    //! Get the complete solution vector inside this object
    /*!
      \param val (output) Vector containing complete solution
    */
    virtual void GetAlgSysVector(CFSVector & val) const = 0;
  
  
    //! Get the pointer to the CFSVector inside this object
    /*!
      \param ptrToVec (output) Pointer to vector inside this object
    */
    virtual void GetAlgSysVectorPointer(CFSVector* &ptrToVec) = 0;


    //! Copies the data from the given pointer

    //! \param ptr (input) Pointer to raw solution data
    //! \note This method is very intrusive and should only be used
    //! when one can ensure that the internal layout of the solution
    //! matches to the one of the given array. This is the case e.g. for
    //! the solution of the algebraic system.
    virtual void CopyFromAlgSysDataPointer(Double * ptr) {
      Error("BaseNodeStoreSol::CopyFromDataPointer() not implemented here",
            __FILE__, __LINE__);
    }
  
    //@{
    //! Set data pointer for the complete solution
 
    //! Set data pointer for the complete solution.
    //! The old data, where the pointer was pointing 
    //! before, will not be deleted!
    /*!
      \param size (input) Size of data
      \param ptr (input) Pointer to raw solution data
    */
    //! \note: This method is very intrusive and should only be used
    //! when one can ensure, that the internal layout of the solution
    //! matches to the one of the given array. This is the case e.g. for
    //! the solution of the algebraic system.
    virtual void SetAlgSysDataPointer( UInt size, Double * ptr ) 
    {Error("BaseNodeStoreSol::SetDataPointer() not implemented here", __FILE__,
           __LINE__);}

    virtual void SetAlgSysDataPointer( UInt size, Complex * ptr ) 
    {Error("BaseNodeStoreSol::SetDataPointer() not implemented here", __FILE__,
           __LINE__);}
    //@}

    //! Get pointer to data in double* - format

    //! Get pointer to data in double* - format.
    //! For entries of type Double, this method simply returns the 
    //! pointer to the internal data.
    //! For entries of type Complex it will convert the data
    //! to a double* array and return the pointer.
    //! The converted array will contain first all real-valued parts of the
    //! the solution and then the complex part.
    virtual Double* GetAlgSysDoublePointer() = 0;


    //! Output function
    virtual void Print(std::ostream& str) = 0;

    /** debug output for developers */
    void Dump()
    {
  	  std::cout << "BaseNodeStoreSol: solutions=" << GetNumSolutions() << " size=" << GetSize() << " nodes=" << GetNumNodes() << " total DOFs=" << GetTotalNumDofs() << std::endl;

      StdVector<SolutionType> solTypes;
      GetSolutionTypes(solTypes);
      
      for(UInt i = 0; i < solTypes.GetSize(); i++) {
          SolutionType type = solTypes[i];
          std::cout << i << ": type=" << type << " dof=" << GetDof(type) << std::endl;
      }  
    };
    


    /////////////////////////////////////////
    // Transformation routines for mapping //
    // pde to mesh solution and vice versa //
    /////////////////////////////////////////

    virtual void GetElemSolution(CFSVector & elemSol,
                                 const EntityIterator& it) const = 0;

    //! 
    virtual void GetElemSolutionAsMatrix(CFSMatrix & elemSol,
                                         const EntityIterator& it) const = 0;
  
    //! maps the local node solution to the coupling nodes
    virtual
    void NodeSolutionToCoupling(CFSVector & couplingSol,
                                const StdVector<UInt>& nodeNumbers) const = 0;

    // ==========================================================
    // DECLARATION OF INTERFACES FOR NON-DOUBLE STORESOL-CLASSES
    // ==========================================================
#ifndef DOXYGEN_SKIP_THIS
  
#define DEFINE_BASENODESTORESOL_FCT(TYPE)                                               \
  virtual void Init(const TYPE val)                                             \
  {Error("BaseNodeStoreSol::Init() not implemented here", __FILE__, __LINE__);} \
  virtual void Get(const UInt nodeNr,                                        \
                   TYPE & ret) const                                            \
  {Error("BaseNodeStoreSol::Get() not implemented here", __FILE__, __LINE__);}  \
  virtual void Get(const SolutionType type,                                     \
                   const UInt nodeNr,                                        \
                   const UInt dof,                                           \
                   TYPE & ret) const                                            \
  {Error("BaseNodeStoreSol::Get not implemented here", __FILE__, __LINE__);}    \
  virtual void Set(const SolutionType type,                                     \
                   const UInt nodeNr,                                        \
                   const UInt dof,                                           \
                   const TYPE val) const                                        \
  {Error("BaseNodeStoreSol::Set not implemented here", __FILE__, __LINE__);}    \
  virtual void Add(const SolutionType type,                                     \
                   const UInt nodeNr,                                        \
                   const UInt dof,                                           \
                   const TYPE val) const                                        \
  {Error("BaseNodeStoreSol::Add() not implemented here", __FILE__, __LINE__);}  \
   
    DEFINE_BASENODESTORESOL_FCT(Complex);

#endif //DOXYGEN_SKIP_THIS 


  protected:
  
    //! Pointer to grid class
    Grid * ptGrid_;

    //! Pointer to equation map
    EqnMap * eqnMap_;

    //! Type of result
    shared_ptr<ResultInfo> result_;   

    //! Element list
    StdVector<shared_ptr<ElemList> >elems_;
  
    //! Number of nodes
    UInt numNodes_;

    //! Number of different solutions
    UInt numSolutions_;
  
    //! Total number of entries

    //! Total number of entries
    //! (= numNodes * sumOf(solDofs_[i]))
    UInt length_;

    //! 
    UInt lengthVector_;
  
    //! Contains order of solution types

    //! Contains order of solution types.
    //! Stores mapping SolutionType <-> position
    std::map<SolutionType,UInt> solTypes_;

    //! Stores offset of soltypes w.r.t. to beginning

    //! Stores the relative position of one
    //! result w.r.t. to the beginning
    //! (depends on the number of dofs
    //! of all previous results)
    std::map<SolutionType,UInt> solOffset_;

    //! Stores offset
    UInt eqnOffset_;

    //! Contains number of dofs for each quantity
    std::map<SolutionType,UInt> solDofs_;
 
    //! Total number of dofs 
    UInt totalDofs_;

    //! Number of equation dofs
    UInt eqnDofs_;

    //! Array for convertin complex data to double

    //! Array for convertin complex data to double
    //! (ref. GetDoublePointer())
    Vector<Double> convertedData_;
  };
  

  // ======================================================
  // INLINE FUNCTIONS
  // ======================================================

  inline UInt BaseNodeStoreSol::GetTotalNumDofs() const
  { 
    ENTER_IFCN( "BaseNodeStoreSol::GetTotalNumDofs" );
    return totalDofs_; 
  }

  inline UInt BaseNodeStoreSol::GetDof(const SolutionType solType) const
  {
    ENTER_IFCN( "BaseNodeStoreSol::GetDof" );
#ifdef CHECK_INDEX
    std::string errMsg;
    if (numSolutions_ > 1 &&
        solType == NO_SOLUTION_TYPE)
      {
        errMsg = "BaseNodeStoresol::GetDof: Solution contains more than ";
        errMsg += "one solutiontype. \n Please specify a solutiontype as ";
        errMsg += "parameter";
        Error(errMsg.c_str(), __FILE__, __LINE__);
      }
#endif
    UInt dof;
    if (numSolutions_ == 1 && solType == NO_SOLUTION_TYPE)
      dof = (*solDofs_.begin()).second;
    else
      dof = (*solDofs_.find(solType)).second;

    std::string warningMsg;
    if (dof == 0)
      {
        warningMsg = "BaseNodeStoresol::GetDof: Either this type ";
        warningMsg += "of solution was not found or it has ";
        warningMsg += "0 degrees of freedom.";
        Error(warningMsg.c_str(), __FILE__, __LINE__);
      }

    return dof;
  }
  
  inline UInt BaseNodeStoreSol::GetNumNodes() const
  {
    ENTER_IFCN( "BaseNodeStoreSol::GetNumNodes" );
    return numNodes_;
  }
  
  inline UInt BaseNodeStoreSol::GetSize() const 
  {
    ENTER_IFCN( "BaseNodeStoreSol::GetSize()" );
    return length_;
  }
  
  inline UInt BaseNodeStoreSol::GetNumSolutions() const 
  {
    ENTER_IFCN( "BaseNodeStoreSol::GetNumSolutions()" );
    return numSolutions_;
  }
/*
  void BaseNodeStoreSol::Dump()
  {
  	  std::cout << "BaseNodeStoreSol: solutions=" << GetNumSolutions() << " size=" << GetSize() << " nodes=" << GetNumNodes() << " total DOFs=" << GetTotalNumDofs() << std::endl;

      StdVector<SolutionType> solTypes;
      GetSolutionTypes(solTypes);
      
      for(UInt i = 0; i < solTypes.GetSize(); i++) {
          SolutionType type = solTypes[i];
          std::cout << i << ": type=" << type << " dof=" << GetDof(type) << std::endl;
      }  
  }*/

} //end of namespace

#endif
