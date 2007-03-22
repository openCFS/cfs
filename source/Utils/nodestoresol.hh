// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NODESTORESOL_2004
#define FILE_NODESTORESOL_2004

#include "basenodestoresol.hh"
#include "Matrix/matrix.hh"
#include "Utils/vector.hh"
#include "Domain/grid.hh"

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
  //! The defining code for the above example ould look like this:
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
  class NodeStoreSol : public BaseNodeStoreSol{
  public:

    //! Default constructor
    NodeStoreSol();


    //! Constructor with given solution layout
    /*!
      \param numNodes (input) Number of nodes
      \param solTypes (input) Vector containing solution types
      \param solDofs (input) Vector containing number of dofs
    */
    //! \note By calling this constructor, one gets an iniaialized
    //! object, ready to use
    NodeStoreSol(const UInt numNodes, 
                 const StdVector<SolutionType> solTypes, 
                 const StdVector<UInt> solDofs);
             

    //! Constructor with given layout for ONE solutiontype
    /*!
      \param numNodes (input) Number of nodes
      \param solTypes (input) Solution type
      \param solDofs (input) Number dofs
    */
    //! \note By calling this constructor, one gets an iniaialized
    //! object, ready to use
    NodeStoreSol(const UInt numNodes,
                 const SolutionType solType,
                 const UInt numDofs);
             

    //! Copy Constructor
    NodeStoreSol(const NodeStoreSol<TYPE> & x);


    //! Destructor
    virtual ~NodeStoreSol();

    //! Hard coded query if values are complex
    virtual bool IsComplex();
  
    //! Deletes all data and layout information

    //! Deletes all data and layout information.
    //! This method should be called  before the layout of the
    //! solution object is modified via SetNumSolution(), SetNumNodes(),
    //! SetNumDofs() or SetSolutionType().
    void Clear();
 
    //! Set Pointer to nodal equation object
    void SetPtrEQNData( EqnMap * eqnMap,
                        Grid * ptGrid );  

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
    void SetNumSolutions(const UInt nSols);
  

    //! Set the number of solution nodes
    /*!
      \param nNodes (input) Number of solution nodes
    */
    //! \note All entries of this object are deleted afterwards
    void SetNumNodes(const UInt nNodes);


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
                         const UInt numSolution = 0);


    //! Set the number of dofs for one solution type

    //! Set the number of dofs for one solution.
    //! If only one type of solution is stored, the second parameter
    //! can be omitted.
    /*!
      \param dof (input) Number of degrees of freedom
      \param solType (input) Type of solution (ref. enum SolutionType)
    */
    //! \note All entries of this object are deleted afterwards
    void SetNumDofs(const UInt dof, 
                    const SolutionType solType = NO_SOLUTION_TYPE);


    //! Get list of stored solution types
    /*!
      \param solTypes (output) Vector of solution types
    */
    void GetSolutionTypes(StdVector<SolutionType> &solTypes) const;


    //! Access operator for nodal result of given dof
    /*!
      \param nodeNr (input) Number of node
      \param dof (input) Dof of nodal result
    */
    //! \note This method may only be called if object contains only
    //! one type of solution.
    TYPE& operator()(const UInt nodeNr, 
                     const UInt dof);


    //! Access operator for nodal result of given dof
    /*!
      \param nodeNr (input) Number of node
      \param dof (input) Dof of nodal result
    */
    //! \note This method may only be called if object contains only
    //! one type of solution.
    TYPE operator()(const UInt node, 
                    const UInt dof) const;


    //! Get vector with one solutiontype for all nodes
    /*!
      \param (input) Solution type (ref. enum SolutionType)
      \param (output) Vector with given solution type)
    */
    void GetGlobalSolVector(const SolutionType solType,
                            CFSVector & val) const;


    //! Set vector with one solution type for all nodes
    /*!
      \param (input) Solution type (ref. enum SolutionType)
      \param (input) Vector with given solution type)
    */
    //! \note The val-vector must match the internal layout
    //! prescribed by the SetDof() and SetSolutionType() methods!
    void SetAlgSysSolVector(const SolutionType solType, 
                            const CFSVector & val) ;


    //! Set all solution types for one node
    /*!
      \param nodeNr (input) Node number for solution
      \param val (input) Vector containing nodal results
    */
    void SetNodalResult(const UInt nodeNr, 
                        const CFSVector &val);

 
    //! Get all solution types for one node
    /*!
      \param nodeNr (input) Node number for solution
      \param val (output) Vector containing nodal results
    */
    void GetNodalResult(const UInt nodeNr, 
                        CFSVector & val) const;


    //! Get vector of one solution type for all nodes of one given dof
    /*!
      \param solType (input) Solution type (ref. enum SolutionType)
      \param dof (input)  Dof of solType
      \param val (output) Vector containing rsults
    */
    void GetGlobalSolVectorSingleDof(const SolutionType solType, 
                                     const UInt dof, 
                                     CFSVector & val) const;
  

    //! Get solution vector for all nodes  of one given dof
    /*!
      \param dof (input)  Dof of solType
      \param val (output) Vector containing rsults
    */
    //! \note This method may only be called if object contains only
    //! one type of solution.
    void GetGlobalSolVectorSingleDof(const UInt dof, 
                                     CFSVector & val) const;


    //! Get single result of given node for given dof
    /*!
      \param nodeNr (input) Node number of result
      \param dof (input) Dof of result
      \param val (output) Result of node nodeNr for given dof
    */
    //! \note This method may only be called if object contains only
    //! one type of solution.
    void Get(const UInt nodeNr, 
             const UInt dof, 
             TYPE & ret) const;

  
    //! Get single result of given solution type, node and dof
    /*! 
      \param solType (input) Solution type (ref. enum SolutionType)
      \param nodeNr (input) Node number of result
      \param dof (input) Dof of result
      \param val (output)  Result of node nodeNr for given dof
    */
    void Get(const SolutionType solType, 
             const UInt nodeNr, 
             const UInt dof, 
             TYPE & ret) const;
  
  
    //! Set a single entry of a given solution type and a given dof
    /*!
      \param solType (input) Solution type (ref. enum SolutionType)
      \param nodeNr (input) Node number of result
      \param dof (input) Dof of result
      \param val (input)  Result of node nodeNr for given solType and dof
    */  
    void Set(const SolutionType solType, 
             const UInt nodeNr, 
             const UInt dof, 
             const TYPE val);


    //! Add value to a single entry of a given solution type and a given dof
    /*!
      \param solType (input) Solution type (ref. enum SolutionType)
      \param nodeNr (input) Node number of result
      \param dof (input) Dof of result
      \param val (input)  Result of node nodeNr for given solType and dof
    */   
    void Add(const SolutionType solType, 
             const UInt nodeNr, 
             const UInt dof, 
             const TYPE val) const;


    ///! Set the complete solution vector inside this object
    /*!
      \param val (input) Vector containing complete solution
    */
    //! \note The layout of the val-vector must match the layout
    //! of the solution prescribed by SetNumSolutions(), SetNumDofs() and
    //! SetSolutionType()!
    void SetAlgSysVector(const CFSVector & val);


    ///! Get the complete solution vector inside this object
    /*!
      \param val (output) Vector containing complete solution
    */
    void GetAlgSysVector(CFSVector & val) const;

  
    //! Get the pointer to the CFSVector inside this object
    /*!
      \param ptrToVec (output) Pointer to vector inside this object
    */
    void GetAlgSysVectorPointer(CFSVector* & ptrToVec);


    //! Get the pointer to the Vector<TYPE> inside this object
    /*!
      \param ptrToVec (output) Pointer to vector inside this object
    */
    void GetAlgSysVectorPointer(Vector<TYPE>* & ptrToVec);


    //! Copies the data from the given pointer
    /*!
      \param ptr (input) Pointer to raw solution data
    */
    //! \note: This method is very intrusive and should only be used
    //! when one can ensure, that the internal layout of the solution
    //! matches to the one of the given array. This is the case e.g. for
    //! the solution of the algebraic system.
    void CopyFromAlgSysDataPointer(Double * ptr);
  
    //! Set data pointer for the complete solution
 
    //! Set data pointer for the complete solution
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
    void SetAlgSysDataPointer( UInt size, TYPE * ptr);
  
    //! Get pointer to data in double* - format

    //! Get pointer to data in double* - format.
    //! For entries of type Double, this method simply returns the 
    //! pointer to the internal data.
    //! For entries of type Complex it will convert the data
    //! to a double* array and return the pointer.
    //! The converted array will contain first all real-valued parts of the
    //! the solution and then the complex part.
    Double * GetAlgSysDoublePointer();


    //! Get the complete solution vector inside this object
    /*!
      \param val (output) Vector containing complete solution
    */
    Vector<TYPE>&  GetAlgSysVector()
    {return data_;}


    /////////////////////////////////////////
    // Transformation routines for mapping //
    // pde to mesh solution and vice versa //
    /////////////////////////////////////////

    void GetElemSolution(CFSVector & elemSol,
                         const EntityIterator& it) const;


    //!
    void GetElemSolutionAsMatrix(CFSMatrix & elemSol, 
                                 const EntityIterator& it) const;


    //! maps the local node solution to the coupling nodes
    void NodeSolutionToCoupling(CFSVector & couplingSol,
                                const StdVector<UInt>& nodeNumbers) const;
                             

    ///////////////
    // Operators //
    ///////////////
  

    //! assignment operator
    NodeStoreSol<TYPE> & operator= (const NodeStoreSol<TYPE> & x);
  
    //! assignment operator for Base-class
    BaseNodeStoreSol & operator= (const BaseNodeStoreSol & x);
  
    //! Standard output operator (of baseclass)
    void  Print(std::ostream& str);
  protected:
  
    //! contains the solution itself
    Vector<TYPE> data_;
    
  };

} //end of namespace

#endif
