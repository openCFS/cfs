// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BASEELEMSTORESOL_2004
#define FILE_BASEELEMSTORESOL_2004

#include <iostream>
#include <map>

#include "Utils/StdVector.hh"
#include "DataInOut/WriteInfo.hh"
#include "General/environment.hh"
#include "Utils/tools.hh"
#include "MatVec/vector.hh"
#include "PDE/eqnMap.hh"


namespace CoupledField{

  // Forward class declarations
  struct Elem;
  class Grid;
  class SingleVector;
  class DenseMatrix;
  template<class TYPE> class Matrix;
  template<class TYPE> class Vector;


  //! This class is the new interface for handling solutions of PDEs

  //! This class is the new interface for handling the 
  //! (multidimensional) solutions of PDEs (instead of 
  //! old Array-class).
  //!
  //! In principle this class is only a wrapper around
  //! a SingleVector class. Additionally it contains
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
  //! mySol.Init() 
  //! \endverbatim
  //! \note An object of the StoreSolution can only used after the Init()-
  //! routine was called, otherwise an error is reported!
  //! \note Although the names of some methods refer to nodes, this class 
  //! also can handle element solutions.
  class BaseElemStoreSol
  {
  public:
  
    //! Default destructor
  
    //! Since this is a base class used for inheritance we give it a virtual
    //! destructor.
    virtual ~BaseElemStoreSol() {};
  
    //! Assignment operator

    //! Assignment operator has to be reimplemented in derived classes
    //! by the help of type casts.
    virtual BaseElemStoreSol & operator= (const BaseElemStoreSol & x) = 0;

    //! Set Pointer to nodal equation object
    virtual void SetPtrEQNData( EqnMap * eqnMap,
                                Grid * ptGrid)
    {EXCEPTION( "Not implemented here" );}
 
    virtual void SetResult( shared_ptr<ResultInfo>& result ) {
      result_ = result; }

    virtual void AddElemList( shared_ptr<ElemList> elems ) {
      elems_.Push_back( elems ); }

    //! Deletes all data and layout information

    //! Deletes all data and layout information.
    //! This method should be called before the layout of the
    //! solution object is modified via SetNumSolution(), SetNumElems(),
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
      {EXCEPTION("BaseElemStoreSol::Init() not implemented here");}  
  
    //! Set the number of different solution types
    /*! 
      \param nSols (input) Number of different solution types
    */
    //! \note All entries of this object are deleted afterwards
    virtual void SetNumSolutions(const UInt nSols) = 0;
  
    //! Set the number of solution elems
    /*!
      \param nElems (input) Number of solution elems
    */
    //! \note All entries of this object are deleted afterwards
    virtual void SetNumElems(const UInt nElems) = 0;

  
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

  
    //! Get number of elems
    inline UInt GetNumElems() const;
  

    //! Get total length of solution vector

    //! Get total length of solution vector (= totalDofs_ * numElems_)
    inline UInt GetSize() const;

  
    //! Get number of solution types
    inline UInt GetNumSolutions() const;

  
    //! Get vectot with one type of solution. The length
    //! of the vector will be (numer of elems in mesh *
    //! degree of freedoms)
    //! 
    /*!
      \param (input) Solution type (ref. enum SolutionType)
      \param (output) Vector with given solution type)
    */
    virtual void GetGlobalSolVector(const SolutionType solType, 
                                    SingleVector & val) const = 0;

  
    //! Set all solution types for one elem
    /*!
      \param elemNr (input) elem number for solution
      \param val (input) Vector containing nodal results
    */
    virtual void SetElemResult(const UInt elemNr,
                               const SingleVector &val) = 0;

  
    //! Get all solution types for one elem
    /*!
      \param elemNr (input) (Global) Elem number for solution
      \param val (output) Vector containing nodal results
    */
    virtual void GetElemResult(const UInt elemNr,
                               SingleVector & val) const = 0;
  

    //! Get vector of one solution type for all elems of one given dof
    /*!
      \param solType (input) Solution type (ref. enum SolutionType)
      \param dof (input)  Dof of solType
      \param val (output) Vector containing rsults
    */
    virtual void GetGlobalSolVectorSingleDof(const SolutionType solType,
                                             const UInt dof,
                                             SingleVector & val) const = 0;
  
  
    //! Get solution vector for all elems of one given dof
    /*!
      \param dof (input)  Dof of solType
      \param val (output) Vector containing rsults
    */
    //! \note This method may only be called if object contains only
    //! one type of solution.
    virtual void GetGlobalSolVectorSingleDof(const UInt dof,
                                             SingleVector & val) const = 0;


    //! Get single result of given elem for given dof
    /*!
      \param elemNr (input) (Global) Elem number of result
      \param dof (input) Dof of result
      \param val (output) Result of elem elemNr for given dof
    */
    //! \note This method may only be called if object contains only
    //! one type of solution.
    virtual void Get(const UInt elemNr, 
                     const UInt dof,
                     Double & val) const
    {EXCEPTION("BaseElemStoreSol::Get() not implemented here");} 

  
    //! Get single result of given solution type, elem and dof
    /*! 
      \param solType (input) Solution type (ref. enum SolutionType)
      \param elemNr (input) (Global) Elem number of result
      \param dof (input) Dof of result
      \param val (output)  Result of elem elemNr for given dof
    */
    virtual void Get(const SolutionType solType, 
                     const UInt elemNr,
                     const UInt dof, 
                     Double & val) const
    {EXCEPTION("BaseElemStoreSol::Get() not implemented here");}

  
    //! Set a single entry of a given solution type and a given dof
    /*!
      \param solType (input) Solution type (ref. enum SolutionType)
      \param elemNr (input) (Global) Elem number of result
      \param dof (input) Dof of result
      \param val (input)  Result of elem elemNr for given solType and dof
    */  
    virtual void Set(const SolutionType solType, 
                     const UInt elemNr,
                     const UInt dof, 
                     const Double val)
    {EXCEPTION("BaseElemStoreSol::Set() not implemented here");}
  

    //! Add value to a single entry of a given solution type and a given dof
    /*!
      \param solType (input) Solution type (ref. enum SolutionType)
      \param elemNr (input) (Global) Elem number of result
      \param dof (input) Dof of result
      \param val (input)  Result of elem elemNr for given solType and dof
    */   
    virtual void Add(const SolutionType solType, 
                     const UInt elemNr,
                     const UInt dof, 
                     const Double val) const
    {EXCEPTION("BaseElemStoreSol::Add() not implemented here");} 

    //! Output function
    virtual void Print(std::ostream& str) = 0;

    /** debug output for developers */
    void Dump()
    {
      std::cout << "BaseElemStoreSol: solutions=" << GetNumSolutions() << " size=" << GetSize() << " total DOFs=" << GetTotalNumDofs() << std::endl;

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

    //!
    virtual void TransformElemSolution(SingleVector & transformedSolution,
                                       Grid * ptGrid,
                                       StdVector<UInt> & mapping) const
    {EXCEPTION("Not implemented here");}
  
  
    //! maps the local element solution to the coupling elements
    virtual void ElemSolutionToCoupling(SingleVector & couplingSol,
                                        const StdVector<Elem*>& elements,
                                        const SingleVector & elemSol) const
    {EXCEPTION("Not implemented here");}  
    

    // ==========================================================
    // DECLARATION OF INTERFACES FOR NON-DOUBLE STORESOL-CLASSES
    // ==========================================================
#ifndef DOXYGEN_SKIP_THIS
  
#define DEFINE_BASEELEMSTORESOL_FCT(TYPE)                                               \
  virtual void Init(const TYPE val)                                             \
  {EXCEPTION("BaseElemStoreSol::Init() not implemented here");} \
  virtual void Get(const UInt elemNr,                                        \
                   TYPE & ret) const                                            \
  {EXCEPTION("BaseElemStoreSol::Get() not implemented here");}  \
  virtual void Get(const SolutionType type,                                     \
                   const UInt elemNr,                                        \
                   const UInt dof,                                           \
                   TYPE & ret) const                                            \
  {EXCEPTION("BaseElemStoreSol::Get not implemented here");}    \
  virtual void Set(const SolutionType type,                                     \
                   const UInt elemNr,                                        \
                   const UInt dof,                                           \
                   const TYPE val) const                                        \
  {EXCEPTION("BaseElemStoreSol::Set not implemented here");}    \
  virtual void Add(const SolutionType type,                                     \
                   const UInt elemNr,                                        \
                   const UInt dof,                                           \
                   const TYPE val) const                                        \
  {EXCEPTION("BaseElemStoreSol::Add() not implemented here");}
    
    DEFINE_BASEELEMSTORESOL_FCT(Complex);

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

    //! Number of elements
    UInt numElems_;

    //! Number of different solutions
    UInt numSolutions_;
  
    //! Total number of entries

    //! Total number of entries
    //! (= numElems * sumOf(solDofs_[i]))
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

    //! Array for convertin complex data to double

    //! Array for convertin complex data to double
    //! (ref. GetDoublePointer())
    Vector<Double> convertedData_;
  };
  

  // ======================================================
  // INLINE FUNCTIONS
  // ======================================================

  inline UInt BaseElemStoreSol::GetTotalNumDofs() const
  { 
    return totalDofs_; 
  }

  inline UInt BaseElemStoreSol::GetDof(const SolutionType solType) const
  {
    if (numSolutions_ == 1 && solType == NO_SOLUTION_TYPE)
      return (*solDofs_.begin()).second;
    else
      return (*solDofs_.find(solType)).second;
  }
  
  inline UInt BaseElemStoreSol::GetNumElems() const
  {
    return numElems_;
  }
  
  inline UInt BaseElemStoreSol::GetSize() const 
  {
    return length_;
  }
  
  inline UInt BaseElemStoreSol::GetNumSolutions() const 
  {
    return numSolutions_;
  }


} //end of namespace

#endif
