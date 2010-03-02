// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_MIXED_EQNMAP_HH
#define FILE_CFS_MIXED_EQNMAP_HH

#include "eqnMap.hh"
#include "disContEqnMap.hh"
//#include "Utils/StdVector.hh"
//#include "Domain/resultInfo.hh"
//#include "Domain/bcs.hh"
//#include "Utils/vector.hh"
//
namespace CoupledField {
  
  //class EntityList;
  //class NodeList;
  //class ElemList;
  //class EntityIterator;
  //class ResultInfo;
  //class InfoNode; 

  //! Class for mapping entities and continuous ansatz functions to equation numbers
  class MixedEqnMap : public EqnMap {
    
  public:
    
    //! Constructor
    MixedEqnMap( Grid* ptGrid, PdeIdType, bool usePenalty);

    //! Destructor
    virtual ~MixedEqnMap();
    
    // ======================================================================
    // SET AND INITIALIZATION METHODS
    // ======================================================================
    
    //@{
    //! \name Set And Initialization Methods
    
    //! Add a new result
    virtual void AddResult( ResultInfo& result, 
                    shared_ptr<EntityList> list );
    
    //! Set the homogeneous boundary conditions
    virtual void SetHomoDirichletBCs( HdBcList& hdBcs );
        
    //! Set the in-homogeneous boundary conditions
    virtual void SetInhomDirichletBCs( IdBcList& idBcs );

    //! Set the constraint conditions
    virtual void SetConstraints( ConstraintList& constraints );
    
    //! Finalize setup and calculate mapping
    virtual void Finalize();

    //! Maps the equation numbers according to the reordering

    //! If the algebraic system has re-ordered the equations, e.g. to improve
    //! the performance of a sparse direct solver, this method must be used to
    //! incorporate the resulting change into the dof to eqn map.
    //! \param order permutation vector containing the re-ordering, i.e.
    //!              order[k-1] is the new equation number for the dof that
    //!              was previously associated to the k-th equation number.
    //! \note
    //! - It is safe to pass a NULL pointer to the method. In this case
    //!   the method assumes that no re-ordering was performed by the
    //!   algebraic system.
    //! - Since the re-ordering is private property of the equation data
    //!   object, we re-set order to NULL in this method.
    //! - For the sake of a low memory foot-print we delete the permutation
    //!   buffer once the reordering was incorporated
    virtual void ReorderMapping( const StdVector<UInt>& order ){
      this->contMap->ReorderMapping(order);
      this->disContMap->ReorderMapping(order);
    }
    
    //@}
    // ======================================================================
    // GENERAL QUERY METHODS
    // ======================================================================

    //@{
    //! \name General Query Methods
    //! Return true if eqns are assigned
    virtual inline bool IsFinalized() const {return (contMap->IsFinalized() && disContMap->IsFinalized());}

    //! Return the total number of equations
    virtual inline UInt GetNumEqns() const { 
      //return contMap->GetNumEqns()+disContMap->GetNumEqns(); 
      return disContMap->GetNumEqns();
      }

    //! Return the equation number of the last unfixed degree of freedom

    //! This method returns the equation number of the last free degree of
    //! freedom. We consider a degree of freedom to be free, if it is not
    //! fixed by either
    //! - homogeneous Dirichlet boundary conditions
    //! - inhomogeneous Dirichlet boundary conditions
    //! - constraints (master - slave dof relations)
    //!
    //! The return value actually depends on the status of the usePenalty_ flag.
    //! If it is true,
    //! then inhomogeneous Dirichlet boundary conditions are treated by the
    //! penalty method and must be considered free dofs as well. Thus, in this
    //! case we return the total number of equations. Otherwise we return only
    //! the number of equations withou a dirichlet boundary condition
    virtual UInt GetNumLastFreeDof() const { 
      UInt eqns = 0;
      eqns = contMap->GetNumLastFreeDof();
      eqns += disContMap->GetNumLastFreeDof()-contMap->GetNumEqns();
      return eqns; 
    }
                              
    //! Return number of real inhomogeneous Dirichlet boundary conditions

    //! This method returns the number of 'real' boundary conditions, i.e.
    //! douplets of nodes are already removed, so this number represents
    //! the number of equations, which are not 'free'
    virtual inline UInt GetNumInHomDirichletEqns() const {
      return contMap->GetNumInHomDirichletEqns()+disContMap->GetNumInHomDirichletEqns();;
    }
    //! Return number of constraint slave equations.

    //! This method returns the number of equations, which are fixed
    //! by a slave constraint condition.
    virtual inline UInt GetNumConstraintSlaveEqns() const {
      return contMap->GetNumConstraintSlaveEqns()+disContMap->GetNumConstraintSlaveEqns();;
    }
                              
    //@}

    // ======================================================================
    // EQUATION MAPPING
    // ======================================================================
    
    //@{
    //! \name Equation Mapping
    //! Map result and entities to equation numbers
    virtual void GetEqns( StdVector<Integer>& eqns,
                  const ResultInfo& result, const EntityIterator& it ) const;

    //! Map result, entities and dof numbers to equation numbers
    virtual void GetEqns( StdVector<Integer>& eqns,
                  const ResultInfo& result, const EntityIterator& it,
                  UInt dof ) const;
    
    //! Mpa combination of result, entity and dof to single equation number
    virtual Integer GetEqn( const ResultInfo& result, 
                    const EntityIterator& it, UInt dof ) const;

    //! Name mapping function for obtaining a nodal equation
    virtual Integer GetNodeEqn( const ResultInfo& result, UInt nodeNr, UInt dof );

    //! Name mapping function for obtaining a nodal equation
   virtual  Integer GetNodeEqn( UInt nodeNr, UInt dof ){
      EXCEPTION("The Mixed Equation Map does not support this fuction");
    }

    virtual void GetNodeEqn( const StdVector<UInt>& nodeNrs, 
                     StdVector<Integer>& eqnNrs ){
      EXCEPTION("The Mixed Equation Map does not support this fuction");
    }


    ////! Get all equation Numbers for a given reult
    //virtual void GetResEqns(  StdVector<Integer>& eqns, const ResultInfo& result ) const;

    //@}
    
    // ======================================================================
    // LOCAL/GLOBAL MAPPING OF MESH ENTITIES
    // ======================================================================
    
    //@{ 
    //! \name Local/Global Mapping of Mesh Entities

    //! Get number of local nodes
    virtual inline UInt GetNumLocalNodes() const {return contMap->GetNumLocalNodes()+disContMap->GetNumLocalNodes();}
    
    //! Get number of local elements
    virtual inline UInt GetNumLocalElems() const {return contMap->GetNumLocalElems()+disContMap->GetNumLocalElems();}

    //! Get number of local edges
    virtual inline UInt GetNumLocalEdges() const {return contMap->GetNumLocalEdges()+disContMap->GetNumLocalEdges();}

    //! Get number of local faces
    virtual inline UInt GetNumLocalFaces() const {return contMap->GetNumLocalFaces()+disContMap->GetNumLocalFaces();}

    //! Map global to local node numbers
    //! (needed for nodal displacement of grid)
    virtual void Mesh2PdeNode(StdVector<UInt> & PdeNodes,
                      const StdVector<UInt> & MeshNodes) const { EXCEPTION(" Mesh2PdeNode is not supported by the mixed Equation map!");}

    //! Map global to local node number
    //! (needed for nodal displacement of grid)
    virtual UInt Mesh2PdeNode(const UInt meshNode) const { EXCEPTION(" Mesh2PdeNode is not supported by the mixed Equation map!");}

    
    //! Map local to global node number
    virtual void Pde2MeshNode(StdVector<UInt> & meshNodes,
                      const StdVector<UInt> & pdeNodes) const { EXCEPTION(" Mesh2PdeNode is not supported by the mixed Equation map!");}

    //! Map local to global node number
    virtual UInt Pde2MeshNode(const UInt pdeNode) const { EXCEPTION(" Mesh2PdeNode is not supported by the mixed Equation map!");}

    //! Map global to local elem number
    virtual UInt Mesh2PdeElem(const UInt elemNumGlob) const { EXCEPTION(" Mesh2PdeNode is not supported by the mixed Equation map!");}

    //! Map local to global elem number
    virtual UInt Pde2MeshElem(const UInt elemNumLoc) const { EXCEPTION(" Mesh2PdeNode is not supported by the mixed Equation map!");}

    //@}

    /** Give all details of the mapping to the info.xml file it triggered 
     * such in the comman line. One can gain similar data via the loggin 
     * (in the debug version) */ 
    virtual void ToInfo(InfoNode* in) const;
    
  private:

    //! Pointer to the equationMaps
    EqnMap* contMap;
    DiscontinuousEqnMap* disContMap;
  };
  
  
#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class EqnMap
  //! 
  //! \purpose
  //!
  //! This class is the base class which defines the interface of all
  //! objects that deal with the assignment of equation numbers to degrees
  //! of freedom. The following gives an overview on the meaning of the
  //! different equation numbers
  //! <center><img src="../AddDoc/splitting.png"></center>
  //!
  //! \collab 
  //!
  //! \implement 
  //!
  //! \status In use
  //!
  //! \unused 
  //!
  //! \improve
  //! 

#endif
}

#endif
