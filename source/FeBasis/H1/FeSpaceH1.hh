#ifndef FILE_CFS_FESPACEH1_HH
#define FILE_CFS_FESPACEH1_HH

#include "FeBasis/FeSpace.hh"

///////////////////////////////////////////////////////////////////
// H1 - Space 
///////////////////////////////////////////////////////////////////

namespace CoupledField {

//!  Class Symbolizing the H1 conforming Finite Element Space (FeSpaceH1) 
/*!
  The Equation Mapping is continuous and lagrange elements are supported
  the order of the Lagrange functions is arbitrary. On the one hand, the standard Linear
  and serendepity (8-Node) quadratic elements, on the other Hand
  Spectral Elements with arbitrary order (will be) are supported.
  The definition of the location of unknowns of the first comes from the 
  preprocessor. For the latter, we create the additional unknowns (and 
  their locations) dynamically. 
  In the current concept, Boundary Conditions for Spectral Approximation
  have to be specified on (Surface) Elements

  TODO: 
    - Adding the additional unknown location to the 
      element shapes to enable correct applicaiton of boundary conditions
      and better postprocessing

*/
class FeSpaceH1 : public FeSpace {


public:

  //! struct for the single equationMap i.e. Nodal Eqns
  //! for the BC map we code the following
  //! <ul> 
  //!   <li>No BC</li>
  //!   <li>HomDirichlet</li>
  //!   <li>InHomDirichlet</li>
  //!   <li>Constraint</li>
  //! </ul>
  struct SingleEqnMap{
    std::map<Integer, StdVector<Integer> > eqns;
    std::map< Integer,StdVector<BcType> > BcKeys;
    
    //! Offset for every 
    //special treatment of constraints we store for every
    //slave node, DOF pair the master node,Dof pair
    //this is very circuitous we need to find an easier way 
    std::map<std::pair<Integer,Integer>, std::pair<Integer,Integer>  > constraintNodes;
    
    StdVector<Integer>& operator[](Integer key){
      return eqns[key];
    }
  };

  //! struct for the single equationMap i.e. Edge,Face etc Eqns
  //! for the BC map we code the following
  //! <ul> 
  //!   <li>No BC</li>
  //!   <li>HomDirichlet</li>
  //!   <li>InHomDirichlet</li>
  //!   <li>Contrain</li>
  //! </ul>
  struct MultiEqnMap{
    
    //! The matrix is arranged as (components/dofs x numFncs)
    std::map< Integer, Matrix<Integer> > eqns;
    
    //! The first vector adresses the components, the inner one the number of functions
    std::map< Integer, StdVector< StdVector<BcType> > > BcKeys;
    //special treatment of constraints we store for every
    //slave DOF the master Dof
    std::map< Integer,  Integer> slaveMasterEqns;
    
    Matrix<Integer>& operator[](Integer key){
      return eqns[key];
    }
  };

  //! Constructor
  FeSpaceH1(PtrParamNode aNode, PtrParamNode infoNode);

  //! Destructor
  ~FeSpaceH1();

  //! Returns the number of (vectorial) unknowns on the element
  virtual UInt GetNumFunctions( const EntityIterator ent );
  
  //! Return equation numbers
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ); 

  //! Return equation numbers for a specific dof
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                        , UInt dof ); 
  //! Return equation numbers for a specific dof and entitytype
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                        , UInt dof, BaseFE::EntityType );
  
  //! get a nodal equation number
  //TODO: in the context of higher order methods and discontiunous methods
  //this should not be used anymore
  //here only for compatibility to nodestoresolution class
  virtual UInt GetNodeEqn(UInt nodeNr, UInt dof);


  //! Get Equation numbers for a specific element
  virtual void GetElemEqns(StdVector<Integer>& eqns,const Elem* elem);

  //! Get Equation numbers for a specific element
  virtual void GetElemEqns(StdVector<Integer>& eqns,const Elem* elem, UInt dof);

  
  //! \copydoc FeSpace::GetOlasMappings
  virtual void GetOlasMappings( StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
                                std::map<UInt,StdVector<std::set<Integer> > >&
                                minorBlocks );
  
  //! Add result
  virtual void AddFeFunction( shared_ptr<BaseFeFunction> fct );


    //! Dump information to screen
  virtual void PrintEqnMap();


protected:

  //! Map Nodal BC Equation NUmbers
  virtual void MapNodalBCs();
  
  //! Map Nodal Equation Numbers
  virtual void MapNodalEqns(UInt phase);

  // ====================================================================
  // PROCESS USER INPUT
  // ====================================================================

  //! Nodal Equation Map
  SingleEqnMap nodeMap_;


};

} // end of namespace
#endif
