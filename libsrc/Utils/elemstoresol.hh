#ifndef FILE_ELEMSTORESOL_2004
#define FILE_ELEMSTORESOL_2004

#include "baseelemstoresol.hh"
#include <Matrix/matrix.hh>
#include <Utils/vector.hh>

namespace CoupledField{

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
  template<class TYPE>
  class ElemStoreSol : public BaseElemStoreSol{
  public:

    //! Default constructor
    ElemStoreSol();


    //! Constructor with given solution layout
    /*!
      \param numElems (input) Number of elems
      \param solTypes (input) Vector containing solution types
      \param solDofs (input) Vector containing number of dofs
    */
    //! \note By calling this constructor, one gets an iniaialized
    //! object, ready to use
    ElemStoreSol(const Integer numElems, 
		 const StdVector<SolutionType> solTypes, 
		 const StdVector<Integer> solDofs);
		 

    //! Constructor with given layout for ONE solutiontype
    /*!
      \param numElems (input) Number of elems
      \param solTypes (input) Solution type
      \param solDofs (input) Number dofs
    */
    //! \note By calling this constructor, one gets an iniaialized
    //! object, ready to use
    ElemStoreSol(const Integer numElems,
		 const SolutionType solType,
		 const Integer numDofs);
    

    //! Copy Constructor
    ElemStoreSol(const ElemStoreSol & x);


    //! Destructor
    ~ElemStoreSol();

    //! Set Pointer to nodal equation object
    void SetPtrEQNData(NodeEQN * ptNodeEQN,
		       Grid * ptGrid,
		       Integer level);
  
    //! Deletes all data and layout information

    //! Deletes all data and layout information.
    //! This method should be called  before the layout of the
    //! solution object is modified via SetNumSolution(), SetNumElems(),
    //! SetNumDofs() or SetSolutionType().
    void Clear();

    //! Initialization of the StoreSolution-object with 0-element(REQUIRED)
  
    //! Initializes the object AFTER the other layout
    //! functions have been called
    //! \note Only after calling Init, the object can
    //! store information
    //! \note Only necessary when default constructor was used
    void Init();

    //! Initialization of the StoreSolution-object (REQUIRED)
  
    //! Initializes the object AFTER the other layout
    //! functions have been called
    //! \note Only after calling Init, the object can
    //! store information
    /*!
      \param val (input) Value the object gets initialized with
    */
    //! \note Only necessary when default constructor was used
    void Init(const TYPE val);


    //! Set the number of different solution types
    /*! 
      \param nSols (input) Number of different solution types
    */
    //! \note All entries of this object are deleted afterwards
    void SetNumSolutions(const Integer nSols);
  

    //! Set the number of solution elems
    /*!
      \param nElems (input) Number of solution elems
    */
    //! \note All entries of this object are deleted afterwards
    void SetNumElems(const Integer nElems);


    //! Set the type of solution and its position 

    //! Set the tpye of solution and its order of occurence
    //! in the solution vector (e.g. in Piezo, the mechanics displacement
    //! is stored previous to the electric potential)
    /*!
      \param solType (input) Type of solution (ref. enum SolutionType)
      \param numSolution (input) Occurence order in the solution vector
    */
    //! \note All entries of this object are deleted afterwards
    void SetSolutionType(const SolutionType solTypes, 
			 const Integer numSolution = 0);


    //! Set the number of dofs for one solution type

    //! Set the number of dofs for one solution.
    //! If only one type of solution is stored, the second parameter
    //! can be omitted.
    /*!
      \param dof (input) Number of degrees of freedom
      \param solType (input) Type of solution (ref. enum SolutionType)
    */
    //! \note All entries of this object are deleted afterwards
    void SetNumDofs(const Integer dof, 
		    const SolutionType solType = NO_SOLUTION_TYPE);


    //! Get list of stored solution types
    /*!
      \param solTypes (output) Vector of solution types
    */
    void GetSolutionTypes(StdVector<SolutionType> &solTypes) const;


    //! Access operator for element result of given dof
    /*!
      \param elemNr (input) Global elem number
      \param dof (input) Dof of nodal result
    */
    //! \note This method may only be called if object contains only
    //! one type of solution.
    TYPE& operator()(const Integer elemNr, 
		     const Integer dof);


    //! Access operator for element result of given dof
    /*!
      \param elemNr (input) (Global) Elem number
      \param dof (input) Dof of nodal result
    */
    //! \note This method may only be called if object contains only
    //! one type of solution.
    TYPE operator()(const Integer elemNr, 
		    const Integer dof) const;

    //! Get vectot with one type of solution. The length
    //! of the vector will be (numer of elems in mesh *
    //! degree of freedoms)
    //! 
    /*!
      \param (input) Solution type (ref. enum SolutionType)
      \param (output) Vector with given solution type)
    */
    void GetGlobalSolVector(const SolutionType solType, 
			    CFSVector & val) const;

  
    //! Set all solution types for one elem
    /*!
      \param elemNr (input) elem number for solution
      \param val (input) Vector containing nodal results
    */
    void SetElemResult(const Integer elemNr,
		       const CFSVector &val);
  
  
    //! Get all solution types for one elem
    /*!
      \param elemNr (input) (Global) Elem number for solution
      \param val (output) Vector containing nodal results
    */
    void GetElemResult(const Integer elemNr,
		       CFSVector & val) const;
  
  
    //! Get vector of one solution type for all elems of one given dof
    /*!
      \param solType (input) Solution type (ref. enum SolutionType)
      \param dof (input)  Dof of solType
      \param val (output) Vector containing rsults
    */
    void GetGlobalSolVectorSingleDof(const SolutionType solType,
				     const Integer dof,
				     CFSVector & val) const;
  
  
    //! Get solution vector for all elems of one given dof
    /*!
      \param dof (input)  Dof of solType
      \param val (output) Vector containing rsults
    */
    //! \note This method may only be called if object contains only
    //! one type of solution.
    void GetGlobalSolVectorSingleDof(const Integer dof,
				     CFSVector & val) const;

    //! Get single result of given elem for given dof
    /*!
      \param elemNr (input) (Global) Elem number of result
      \param dof (input) Dof of result
      \param val (output) Result of elem elemNr for given dof
    */
    //! \note This method may only be called if object contains only
    //! one type of solution.
    void Get(const Integer elemNr, 
	     const Integer dof,
	     TYPE & val) const;
  

    //! Get single result of given solution type, elem and dof
    /*! 
      \param solType (input) Solution type (ref. enum SolutionType)
      \param elemNr (input) (Global) Elem number of result
      \param dof (input) Dof of result
      \param val (output)  Result of elem elemNr for given dof
    */
  
    void Get(const SolutionType solType, 
	     const Integer elemNr,
	     const Integer dof, 
	     TYPE & val) const;


    //! Set a single entry of a given solution type and a given dof
      /*!
	\param solType (input) Solution type (ref. enum SolutionType)
	\param elemNr (input) (Global) Elem number of result
	\param dof (input) Dof of result
	\param val (input)  Result of elem elemNr for given solType and dof
      */  
      void Set(const SolutionType solType, 
	       const Integer elemNr,
	       const Integer dof, 
	       const TYPE val);
 

      //! Add value to a single entry of a given solution type and a given dof
      /*!
	\param solType (input) Solution type (ref. enum SolutionType)
	\param elemNr (input) (Global) Elem number of result
	\param dof (input) Dof of result
	\param val (input)  Result of elem elemNr for given solType and dof
      */   
      void Add(const SolutionType solType, 
	       const Integer elemNr,
	       const Integer dof, 
	       const TYPE val) const;


    /////////////////////////////////////////
    // Transformation routines for mapping //
    // pde to mesh solution and vice versa //
    /////////////////////////////////////////

    //!
    void TransformElemSolution(CFSVector & transformedSolution,
			       Grid * ptGrid,
			       const Integer level) const;
  

    //! maps the local element solution to the coupling Elems
    void ElemSolutionToCoupling(CFSVector & couplingSol,
				const StdVector<Elem*>& elements,
				const CFSVector & elemSol) const;
  

    ///////////////
    // Operators //
    ///////////////
  

    //! assignment operator
    ElemStoreSol & operator= (const ElemStoreSol & x);
  
    //! assignment operator for Base-class
    BaseElemStoreSol & operator= (const BaseElemStoreSol & x);
  
    //! Standard output operator (of baseclass)
    void  Print(std::ostream& str);
  protected:
  
    //! contains the solution itself
    Vector<TYPE> data_;
    
};

 
  // ======================================================
  // EXPLICIT TEMPLATE INSTANTIATION
  // ======================================================
#if defined(__GNUC__) 
  template class ElemStoreSol<Double>;
  template class ElemStoreSol<Complex>;
  //template class ElemStoreSol<Integer>;
#endif 

} //end of namespace

#endif
