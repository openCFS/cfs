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
    virtual Elem::FEType feType() {return feType_;}
    
    //! Get total number of dofs for particular dof
    virtual UInt GetNumFncs( ) {return actNumFcns_;}

    //! Get number of shape functions for a given type
    virtual void GetNumFncs( StdVector<UInt>& numFcns,
                             ElementEntityType fctEntityType,
                             UInt dof = 1) = 0;

    //! Get the permutation Vector for a given Face or Edge
    //! e.g. If asked for a face, the element will check the flags
    //! of this face and return a vector of size NumberOfFncs on the Face
    //! holding the correct ordering 
    /*!
      \param fncPermutation (output) The Permuation Vector 
      \param ptElem (input) pointer to Grid Element to get grip of flags 
      \param fctEntityType (input) The Entity type, Node/Edge/Face
      \param entNumber (input) The local entity number 
    */
    virtual void GetFncPermutation( StdVector<UInt>& fncPermutation,
                                    const Elem* ptElem,
                                    ElementEntityType fctEntityType,
                                    UInt entNumber){
    };

    //! Get value of all shape fnc at integration point ip
    /*!
      \param S (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
      \param ip (input) Integration point
    */
    virtual void GetShFnc( Vector<Double>& S, const LocPoint& lp,
                              const Elem* ptElem,  UInt comp = 1 ){}

    
    //! Get local derivatives of all shape fnc at arbitrary local point
    //! Local means here on the reference element
    /*! 
      \param S (output) Matrix with global derivatives of all shape functions
      \f [ \left( \begin{array}{ccc} N_{1,dx} & N_{1,dy} & \cdots \\
      N_{2,dx} & N_{2,dy} & \cdots \\
      \cdots     & \cdots      & \cdots \end{array}\right) \f ]
      \param LCoord (input) Local Coordinates of evalutaion point
      \param CornerCoords (input) Coordinates of element corners
      \f [ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
      \cdots & \cdots & \cdots \end{array} \right) \f ]       
    */
    virtual  void GetDerivShFnc(Matrix<Double> & Deriv, 
                                    const LocPoint& lp,
                                    const Elem * elem, 
                                    UInt comp = 1){}


    //! Set the isotropic order of the Element. This methods gets overwritten 
    //! by the child classes to calculate the number of functions according to
    //! the given order
    //! \param order (input) The desired order of the element
    virtual void SetIsoOrder(UInt order){};

    //! Set the Anisotropic order of the Element. This methods gets overwritten 
    //! by the child classes to calculate the number of functions according to
    //! the given order
    //! \param order (input) vector of element orders for each space direction 
    virtual void SetAnIsoOrder(StdVector<UInt> order){};

    //! Calculates corresponding volume point of neighbouring surfaces
    //! For a given surface element and a neighbouring volume element this
    //! method calculates the local volume-coordinates out of the given
    //! local surface-coordinates, which have one less dimension.
    //! This can be used to get the corresponding volume coordinates of
    //! the integration points of a surface. Therefore it calculates
    //! on which side of the volume element the surface elemente lies
    //! and creates the according volume point.
    /*!
      \param surfConnect (input) Node numbers of surface element
      \param volConnect (input) Node numbers of colume element
      \param surfIntPoint (input) Surface integration point, which gets mapped
      onto the volume element
      \param volIntPoint (output) Corresponding volume integration point
    */
//    virtual void GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
//                                                const StdVector<UInt> & volConnect,
//                                                 const Vector<Double> & surfIntPoint,
//                                                Vector<Double> & volIntPoint)
//    {
//      EXCEPTION( "BaseFE::GetLocalIntPoints4Surface: Not implemented here" );}


  protected:

     //! Define variables of this class
     virtual void Init()
     { EXCEPTION("Init not implemented for this element type! ");};
     
     
     // =======================================================================
     //  PRE CALCULATION OF SHAPE FUNCTIONS AT INTEGRATION POINTS
     // =======================================================================
     
     //! Set value of shape fnc at integration points
     virtual void SetShapeFncAtIp();


     //! Set value of shape fnc derivatives at integration points
     virtual void SetShapeFncDerivAtIp();    

    //! Actual number of ansatz functions
    UInt actNumFcns_;
    
    //! Actual Order of the element
    UInt order_;

    Elem::FEType feType_;
  };

}

#endif // FILE_BASEFE_2003
