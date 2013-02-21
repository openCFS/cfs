// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_BCS_HH
#define FILE_CFS_BCS_HH

#include <set>
#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

  // forward class declaration
  class EntityList;
  class EqnMap;
  struct ResultInfo;
  template <class TYPE> class StdVector;

  //! Definition of a homogeneous Dirichlet boundary condition
  struct HomDirichletBc {

    //! Constructor
    HomDirichletBc();

    //! Destructor
    virtual ~HomDirichletBc();

    //! EntityList where the condition is defined on
    shared_ptr<EntityList> entities;

    //! Type of result the boundary condition is assigned with
    shared_ptr<ResultInfo> result;

    //! Set of indices, to which the Dirichlet BC applies
    std::set<UInt> dofs;

    //! Output content to info.xml 
    virtual void ToInfo(PtrParamNode in) const;
  };


  // -------------------------------------------------------------------------

  // Inhomogeneous Dirichlet boundary condition
  struct InhomDirichletBc : public HomDirichletBc {

    //! Constructor
    InhomDirichletBc();

    //! Coefficient function for the values
    PtrCoefFct value; 

    //!  Output our content to info.xml 
    virtual void ToInfo(PtrParamNode in) const;
    
    //! Flag, if updated geometry is to be used
    bool updatedGeo;
  };

  // -------------------------------------------------------------------------
  struct Constraint {

    //! Constructor
    Constraint();

    //! Master entitiyList
    shared_ptr<EntityList> masterEntities;

    //! Slave entityList
    shared_ptr<EntityList> slaveEntities;

    //! Degree of freedom for master entities
    UInt masterDof;

    //! Degree of freedom for slave entities
    UInt slaveDof;

    //! Type of result the boundary condition is assigned with
    shared_ptr<ResultInfo> result;

    //! Equation map
    shared_ptr<EqnMap> eqnMap;

    /** does this constraint originate from periodic bcs */
    bool periodic;
  };

  // -------------------------------------------------------------------------

  // Public typedefs
  typedef StdVector<shared_ptr<HomDirichletBc> > HdBcList;
  typedef StdVector<shared_ptr<InhomDirichletBc> > IdBcList;
  typedef StdVector<shared_ptr<Constraint> > ConstraintList;


}
#endif
