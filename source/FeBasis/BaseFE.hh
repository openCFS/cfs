#ifndef FILE_BASEFE_HH
#define FILE_BASEFE_HH

#include <map>

#include "General/Environment.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Utils/StdVector.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Forms/IntScheme.hh"

namespace CoupledField
{
  // Forward class declaration
  class AlgebraicSys;

  // ========================================================================
  //  FINITE ELEMENT CLASS
  // ========================================================================
  
  //! Base class for computational element 
  
  //! This class serves as an interface for all kinds of computational
  //! elements. By computational element, we denote that this
  //! class (and especially the derived ones) provide methods for calculating
  //! the basis / shape functions adn their derivatives. 
  //! Depending on the geometric shape (line, quad, hex ...), the type of 
  //! continuity of the scalar / vector functions (H1, H(curl), H(div), L2) 
  //! and the type of polynomials used (Lagrange, hierarchical Legendre / 
  //! Jacobi), we have derived classes.
  class BaseFE {
  public:
    
    //! Defines available element entity types
    
    //! Define element entity types, on which unknowns are associated with,
    //! i.e. vertex (= corner nodes), node (= all nodes), edges, faces and
    //! the interior.
    typedef enum { ALL, VERTEX, NODE, EDGE, FACE, INTERIOR} EntityType;
    
    //! Static Enum for conversion of ElemShapeType
    static Enum<EntityType> entityType;

    
    //! Constructor (does nothing)
    BaseFE();

    //! Copy constructor (does something)
    BaseFE(const BaseFE & fe);

    //! declare clone operation
    virtual BaseFE* Clone()=0;

    //! Decstructor
    virtual ~BaseFE();
    
    //! Element shape of reference element
    ElemShape shape_;
    
    // ========================================================================
    //  B A S I S   F U N C T I O N S
    // ========================================================================
    
    //! Return FE-Type
    virtual Elem::FEType FeType() {return feType_;}
    
    //! Flag, if element has isotropic polynomial order
    virtual bool IsIsotropic() {
      return true; 
    }
    
    //! Get total number of dofs for particular dof
    virtual UInt GetNumFncs( ) {return actNumFncs_;}
    
    //! Get accumulated number of functions for entities of given type 
    virtual UInt GetNumFncs( EntityType entType,
                             UInt dof = 1);

    //! Get number of functions for all entities of given type (edge/face...)
    
    //! This method delivers for a given entity type (NODE/EDGE/FACE/INNER) 
    //! the number of unknowns associated with each NODE/EDGE/FACE/INNER.
    //! \param numFncs (out) Vector containing the number of unknowns associated
    //!                      with each of given type
    //! \param entType (in) Type of entity (NODE, EDGE, FACE, INNER)
    //! \Note: For Lagrangian elements, only the entity type NODE is defined.
    virtual void GetNumFncs( StdVector<UInt>& numFcns,
                             EntityType entType,
                             UInt dof = 1) = 0;
                             
    //@{ \name Standard bilinear shape function 
    //! Get value of all shape fnc at local poiont lp
    /*!
    \param S (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
    \param ip (input) Integration point
    */
    virtual void GetShFnc( Vector<Double>& S, const LocPoint& lp,
                           const Elem* ptElem,  UInt comp = 1 ) {
     EXCEPTION("This GetShFnc is not implemented for BaseFE")
    }
    
    //! Return global curl of shape functions
    virtual void GetCurlShFnc( Matrix<Double>& curl, 
                               const LocPointMapped& lp,
                               const Elem* elem, UInt comp = 1 ) {
      EXCEPTION("This GetCurlShFnc is not implemented for BaseFE")
    }   


    //! Flag, if element has true nodal permutation
    virtual bool NeedsNodalPermutation() {
      return false;
    }
    
    
    //! Get the nodal permutation Vector for a given Face or Edge
    //! e.g. If asked for a face, the element will check the flags
    //! of this face and return a vector of size NumberOfFncs on the Face
    //! holding the correct ordering 
    /*!
      \param fncPermutation (output) The permutation vector 
      \param ptElem (input) Pointer to Grid Element to get grip of flags 
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

    //! Obtain isotropic  order of the current element

    //! If there is no order is set, we return 0
    virtual UInt GetIsoOrder() const {
      EXCEPTION("Not implemented");
      return 0;
    }

    //! Return the maximum polynomial order of the element
    virtual UInt GetMaxOrder() const {
      EXCEPTION("Not implemented");
      return 0;
    }
    
    //! Return the maximum polynomial order in each local direction of the element
    virtual void GetAnisoOrder(StdVector<UInt>& order ) const {
      // provide default implementation for isotropic elements
      order.Resize( Elem::shapes[feType_].dim);
      order.Init( GetIsoOrder() );
    }

    //! Returns connect indices for a triangulated simplex
    //! \return Vector of triangle/tet connectivities relative to the reference element 0-based
    virtual void Triangulate(StdVector< StdVector<UInt> > & triConnect){
      Exception("Triangulation not available for the selected element type!");
    }

    //! Compute the coefficient and exponent matrix for an alternate desciption 
    //! of higher order results. I.e. monomial representation.
    //! We assume a test function \f$ \phi_j(\xi,\eta,\zeta)\f$ whihc can be 
    //! expressed as a monomial in the following way
    //! \f[
    //!   \phi_j = \displaystyle \sum_{k=1}^{d} p_k \Phi_{ik}
    //! \f]
    //! In which the coefficients \f$p_k\f$ can be computed by the matrix P 
    //!  \f$d\times3\f$ )which denotes the
    //! exponents of the local coordinates. i.e.
    //! \f[
    //!   p_i = \xi^{P_{i1}} \eta^{P_{i2}} \zeta^{P_{i3}} \ .
    //! \f]
    //! The \f$\Phi_{ik}\f$ are the coefficients of the polynomial and are 
    //! stored in the matrix C.
    //! In case of Lagrange Polynomials C can be found by computing the inverse 
    //! of the vandermonde-matrix.
    //! \param P Matrix storing the exponentials of the local directions
    //! \param C Matrix monomial coefficients
    virtual void ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C){
      Exception("Not implemented in base class");
    }

    virtual std::string GetFeSpaceName(){
      Exception("Not implemented in base class");
      std::string r = "";
      return r;
    }

  protected:

    //! Actual number of basis functions associated with this element
    UInt actNumFncs_;

    //! Geometric type of finite element (line, quad, hex ...)
    Elem::FEType feType_;
    

    //! Flag if can pre-compute shape functions

    //! This flag indicated, if the shape functions and local derivatives can
    //! be precomputed, i.e. they are independent of the global numbering.
    //! This is the case for Lagrangian based shape functions, but nor for
    //! hierarchical shape functions, as they depend on the global connectivity
    //! of the element, to ensure continuous shape functions across element 
    //! boundaries
    bool preComputShFnc_;

  };

}

#endif // FILE_BASEFE_2003
