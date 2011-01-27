// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_BCS_HH
#define FILE_CFS_BCS_HH

#include <ostream>
#include "General/environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

  // forward class declaration
  class EntityList;
  class EqnMap;
  class ResultInfo;
  template <class TYPE> class StdVector;

  //! Definition of a homogeneous Dirichlet boundary condition
  struct HomDirichletBc {

    //! Constructor
    HomDirichletBc();

    /** Add this virtual destructor because we hava a virtul Dump now */
    virtual ~HomDirichletBc();

    //! EntityList where the condition is defined on
    shared_ptr<EntityList> entities;

    //! Type of result the boundary condition is assigned with
    shared_ptr<ResultInfo> result;

    //! Equation map
    shared_ptr<EqnMap> eqnMap;

    //! Degree of freedom index
    UInt dof;

    /** Ouptut our content to info.xml */
    virtual void ToInfo(PtrParamNode in) const;

    /** Just a simple Dump() for developers */
    virtual std::string ToString();
  };


  // -------------------------------------------------------------------------

  // Inhomogeneous Dirichlet boundary condition
  struct InhomDirichletBc : public HomDirichletBc {

    //! Constructor
    InhomDirichletBc();

    //! Value of entities
    std::string value;

    //! Phase value of entities
    std::string phase;

    /** Ouptut our content to info.xml */
    virtual void ToInfo(PtrParamNode in) const;

    virtual std::string ToString();
  };

  // -------------------------------------------------------------------------

  // Inhomogeneous Dirichlet boundary condition read from file
  struct InhomDirichFileBc : public HomDirichletBc {

    //! Constructor
    InhomDirichFileBc();

    //! name of file id in which data is stored
    std::string inputId;

    /** Ouptut our content to info.xml */
    virtual void ToInfo(PtrParamNode in) const;

    virtual std::string ToString();
  };

  // -------------------------------------------------------------------------

  // Inhomgogeneous Neumann boundary condition
  struct InhomNeumannBc : public HomDirichletBc {

    //! Constructor
    InhomNeumannBc();

    //! Value of entities
    std::string value;

    //! Phase value of entities
    std::string phase;

  };


  // -------------------------------------------------------------------------

  // Load definition
  struct LoadBc : public HomDirichletBc {

    //! Constructor
    LoadBc();

    //! Value of entities
    std::string value;

    //! Phase value of entities
    std::string phase;

    /** For multiple excitation optimization */
    std::string weight;

    /** Ouptut our content to info.xml */
    virtual void ToInfo(PtrParamNode in) const;

    virtual std::string ToString();
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
  typedef StdVector<shared_ptr<InhomDirichFileBc> > IdFileBcList;
  typedef StdVector<shared_ptr<InhomNeumannBc> > InBcList;
  typedef StdVector<shared_ptr<LoadBc> > LoadList;
  typedef StdVector<shared_ptr<Constraint> > ConstraintList;


}
#endif
