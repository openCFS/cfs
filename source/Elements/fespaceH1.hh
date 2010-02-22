#ifndef FILE_CFS_FESPACEH1_HH
#define FILE_CFS_FESPACEH1_HH

#include "Elements/fespace.hh"

///////////////////////////////////////////////////////////////////
// H1 - Lower Order / non hierarchical
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
    - Support for Spectral Elements
    - Adding the additional unknown location to the 
      element shapes to enable correct applicaiton of boundary conditions
      and better postprocessing

*/
class FeSpaceH1 : public FeSpace {


public:

  //! Type of basis used
  typedef enum {LAGRANGE, DUAL} BasisType; 


  //! struct for the single equationMap i.e. Nodal Eqns
  //! for the BC map we code the following
  //! <ul> 
  //!   <li>No BC</li>
  //!   <li>HomDirichlet</li>
  //!   <li>InHomDirichlet</li>
  //!   <li>Contrain</li>
  //! </ul>
  struct SingleEqnMap{
    std::map<Integer, StdVector<Integer> > eqns;
    std::map< Integer,StdVector<FeSpace::BC_Type> > BcKeys;
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
    std::map< Integer, Matrix<Integer> > eqns;
    //There is an linker error if we use a matrix for this
    //so we use a vector of vectors
    std::map< Integer, StdVector< StdVector<FeSpace::BC_Type> > > BcKeys;
    //special treatment of constraints we store for every
    //slave DOF the master Dof
    std::map< Integer,  Integer> slaveMasterNodes;
    //Optionally store the connectivity information in a map
    //which associates the e.g. edge number to a set of elements
    //std::map< Integer, StdVector<Integer> > elements;
    Matrix<Integer>& operator[](Integer key){
      return eqns[key];
    }
  };

  //! Constructor
  FeSpaceH1();

  //! Destructor
  ~FeSpaceH1();

  //! Set type of basis
  void SetBasis( BasisType type );

  //! Return pointer to reference element
  virtual BaseFE* GetFe( const EntityIterator ent );

  //! Returns the number of (vectorial) unkowns on the element
  virtual UInt GetNumFunctions( const EntityIterator ent );
  
  //! Return equation numbers
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ); 

  //! Return equation numbers for a specific dof
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                        , UInt dof ); 
  
  //! get a nodal equation number
  //TODO: in the context of higher order methods and discontinous methods
  //this should not be used anymore
  //here only for compatibility to nodestoresolution class
  virtual UInt GetNodeEqn(UInt nodeNr, UInt dof);


  //! Get Equation numbers for a specific element
  virtual void GetElemEqns(StdVector<Integer>& eqns,const Elem* elem);

  //! Get Equation numbers for a specific element
  virtual void GetElemEqns(StdVector<Integer>& eqns,const Elem* elem, UInt dof);

  //! Add result
  virtual void AddFeFunction( shared_ptr<BaseFeFunction> fct );

  //! Map equations i.e. intialize object
  virtual void Finalize();

  //! Reorder the equation Map (just for comptibility)
  virtual void ReorderEqnMap( StdVector<UInt> newOrder );

  virtual void PrintEqnMap();
protected:

  //! Map BC Equation NUmbers
  virtual void MapBCs();

  //! Map Nodal Equation Numbers
  virtual void MapNodalEqns(UInt phase);

  //! Map Edge Equation Numbers
  virtual void MapEdgeEquations(UInt phase);

  //! Map Face Equation Numbers
  virtual void MapFaceEquations(UInt phase);

  //! Map Inerior Equation Numbers
  virtual void MapInteriorEquations(UInt phase);

  //! Type of basis
  BasisType basis_;

  //! Nodal Equation Map
  SingleEqnMap nodeMap_;


};

}
#endif
