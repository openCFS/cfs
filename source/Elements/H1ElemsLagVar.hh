// =====================================================================================
// 
//       Filename:  H1LagrangeVar.hh
// 
//    Description: Lagrangian based elements of arbitrary order (tensor Product elements) 
//
//                 This serves as a base class for all Lagrangian bases Elements of 
//                 arbitrary order. Right now, only lines, quads and hex elements can be 
//                 found here.
//                 The shape functions etc. Will be calculated on request during initialization
//                 depending on the order of the FE-Space
// 
// 
//        Version:  0.1
//        Created:  01/20/2010 09:34:37 AM
//       Revision:  none
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef  FILE_CFS_H1LAGRANGEVAR_HH
#define  FILE_CFS_H1LAGRANGEVAR_HH

#include	"H1Elems.hh"

namespace CoupledField {

  class FeH1LagrangeVar : public FeH1 {

  public:
    //! Constructor
    FeH1LagrangeVar();

      //! Destructor
    virtual ~FeH1LagrangeVar();

    //! Pre-calculate values at integration points
    void SetIntPoints( StdVector<LocPoint>& intPoints );

    //! Get number of shape functions for a given type (NODE/EDGE/FACE/ELEM)
    //! WARNING: if at some point other elements like triangles and tetras are
    //!          available, this general implementation is no langer valid!
    //void GetNumFncs( StdVector<UInt>& numFcns,
    //                 const shared_ptr<AnsatzFct>& fcnType,
    //                 AnsatzFct::FctEntityType fctEntityType,
    //                 UInt dof = 1 );

    //! Get number of shape functions for a given type (NODE/EDGE/FACE/ELEM)
    void GetNumFncs( StdVector<UInt>& numFcns,
                     EntityType fctEntityType,
                     UInt dof = 1 );

    //! This holds only for line,quad and hex,
    //! if other types are available the has to be reimplemented!
    //! Get the permutation Vector for a given Face or Edge
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
    virtual void GetNodalPermutation( StdVector<UInt>& fncPermutation,
                                      const Elem* ptElem,
                                      EntityType fctEntityType,
                                      UInt entNumber);

    ////! Return shape function
    //void GetShFnc( Vector<Double> & S, const LocPoint& lp,
    //               const Elem* ptElem,  UInt comp = 1 );

    ////! Return local derivative of shape function
    //void GetDerivShFnc( Matrix<Double> & deriv, 
    //                    const LocPoint& lp,
    //                    const Elem * elem, 
    //                    UInt comp = 1 );

    //! returns the number of functions for a single edge or face
    UInt GetNumFncsPerEntType( EntityType fctEntityType, UInt dof = 1);

  protected:

    //! Compute shape function at given position
    virtual void CalcShFnc( Vector<Double>& shape,
                            const Vector<Double>& point ) = 0;

    //! Compute local derivative at of shape function at given position
    virtual void CalcDerivShFnc( Matrix<Double> & deriv, 
                                 const Vector<Double>& point ) = 0;

    //! Valued of shape functions at integration points
    StdVector<Vector<Double> > shapeAtIp_;

    //! Value of local derivatives of shape functions at integration points
    StdVector<Matrix<Double> > shapeDerivAtIp_;

  };

  //! Lagrangian line element of variable order
  class FeH1LagrangeLineVar : public FeH1LagrangeVar {

    public:

      //! Constructor
      FeH1LagrangeLineVar();

      //! Destructor
      virtual ~FeH1LagrangeLineVar();

      //! Set the isotropic order of the Element. This methods gets overwritten 
      //! by the child classes to calculate the number of functions according to
      //! the given order
      //! \param order (input) The desired order of the element
      virtual void SetIsoOrder(UInt order);

    protected:

      //! @see FeH1LagrangeExpl::CalcShapeFnc 
      void CalcShFnc( Vector<Double>& shape,
                      const Vector<Double>& point );

      //! @see FeH1LagrangeExpl::CalcDerivShFnc 
      void CalcDerivShFnc( Matrix<Double> & deriv, 
                           const Vector<Double>& point );
  };
  //! Lagrangian quadrilateral element of variable order 
  class FeH1LagrangeQuadVar : public FeH1LagrangeVar {

    public:

      //! Constructor
      FeH1LagrangeQuadVar();

      //! Destructor
      virtual ~FeH1LagrangeQuadVar();

      //! Set the isotropic order of the Element. This methods gets overwritten 
      //! by the child classes to calculate the number of functions according to
      //! the given order
      //! \param order (input) The desired order of the element
      virtual void SetIsoOrder(UInt order);

    protected:

      //! @see FeH1LagrangeExpl::CalcShapeFnc 
      void CalcShFnc( Vector<Double>& shape,
                      const Vector<Double>& point );

      //! @see FeH1LagrangeExpl::CalcDerivShFnc 
      void CalcDerivShFnc( Matrix<Double> & deriv, 
                           const Vector<Double>& point );
   };
   //! Lagrangian hexahedral element of varaiable order
  class FeH1LagrangeHexVar : public FeH1LagrangeVar {

    public:

      //! Constructor
      FeH1LagrangeHexVar();

      //! Destructor
      virtual ~FeH1LagrangeHexVar();

      //! Set the isotropic order of the Element. This methods gets overwritten 
      //! by the child classes to calculate the number of functions according to
      //! the given order
      //! \param order (input) The desired order of the element
      virtual void SetIsoOrder(UInt order);
    protected:

      //! @see FeH1LagrangeExpl::CalcShapeFnc 
      void CalcShFnc( Vector<Double>& shape,
                      const Vector<Double>& point );

      //! @see FeH1LagrangeExpl::CalcDerivShFnc 
      void CalcDerivShFnc( Matrix<Double> & deriv, 
                           const Vector<Double>& point );
  };  
}
#endif     // -----  not FILE_CFS_H1LAGRANGEVAR_HH  -----
