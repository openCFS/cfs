// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BASEFE_2003
#define FILE_BASEFE_2003

#include <map>

#include "General/environment.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Utils/StdVector.hh"
//#include "Domain/ansatzFct.hh"
#include "Domain/elem.hh"
#include "Elements/integrationScheme.hh"

namespace CoupledField
{
  // Forward class declaration
  class BaseSystem;

  //! Base class for description of elements

  //! From this class all special finite elements as triangles, quadrilaterals,
  //! tetraheader, etc. are derived. Some methods are purely virtual.
  //! ToDO: 
  //!       - Get rid of all ..AtIp.. things
  //!       - Get rid of all ..Global.. stuff -> performed in corresponding base element
  //! Note: The base class provides methods for delivering the shape functions for 
  //!       scalar, as well as vectorial shape functions
  class BaseFE {

  public:
    // Due to (Re-) initialization of elements (e.g. in reduced integration)
    friend class Grid;
#ifdef INTEGLIB
    friend class ElemIntegr;
#endif

    
    //! Defines available element entity types
    
    //! Define element entity types, on which unknowns are associated with,
    //! i.e. vertex (= corner nodes), node (= all nodes), edges, faces and
    //! the interior.
    typedef enum { ALL, VERTEX, NODE, EDGE, FACE, INTERIOR} EntityType;
    
    //! Static Enum for conversion of ElemShapeType
    static Enum<EntityType> entityType;
    
    //! constructor (does nothing)
    BaseFE();

    //! decstructor
    virtual ~BaseFE();
    
    //! Element shape of reference element
    ElemShape shape_;
    
    // ========================================================================
    //  B A S I S   F U N C T I O N S
    // ========================================================================
    
    //! Return FE-Type
    virtual Elem::FEType FeType() {return feType_;}
    
    //! Get total number of dofs for particular dof
    virtual UInt GetNumFncs( ) {return actNumFncs_;}
    
    //! Get cumulated number of shape functions for entities of given type (node/edge/face/inner)
    virtual UInt GetNumFncs( EntityType entType,
                             UInt dof = 1);

    //! Get number of shape functions for all entities of given type (edge/face etc.)
    
    //! This method delivers for a given entity type (NODE/EDGE/FACE/INNER) the number
    //! of unknowns associated with each NODE/EDGE/FACE/INNER.
    //! \param numFncs (out) Vector containing the number of unkowns associated with each 
    //!                      of given type
    //! \param entType (in) Type of entity (NODE, EDGE, FACE, INNER)
    //! \Note: For Lagrangian elements, only the entity type NODE is defined.
    virtual void GetNumFncs( StdVector<UInt>& numFcns,
                             EntityType entType,
                             UInt dof = 1) = 0;

    //! Get the nodal permutation Vector for a given Face or Edge
    //! e.g. If asked for a face, the element will check the flags
    //! of this face and return a vector of size NumberOfFncs on the Face
    //! holding the correct ordering 
    /*!
      \param fncPermutation (output) The Permuation Vector 
      \param ptElem (input) pointer to Grid Element to get grip of flags 
      \param fctEntityType (input) The Entity type, Node/Edge/Face where the 
                                   nodes are located at
      \param entNumber (input) The local entity number 
    */
    //! \Note: For non-Lagrangian elements, only the VERTEX fctEntityType returns
    //!        values, as the higher order coefficients are not associated to nodes.
    virtual void GetNodalPermutation( StdVector<UInt>& fncPermutation,
                                      const Elem* ptElem,
                                      EntityType fctEntityType,
                                      UInt entNumber){
    };

    //! Set the isotropic order of the Element. This methods gets overwritten 
    //! by the child classes to calculate the number of functions according to
    //! the given order
    //! \param order (input) The desired order of the element
    virtual void SetIsoOrder(UInt order){};

  protected:

    //! Actual number of ansatz functions
    UInt actNumFncs_;

    //! Geometrix type of finite element
    Elem::FEType feType_;
  };

}

#endif // FILE_BASEFE_2003
