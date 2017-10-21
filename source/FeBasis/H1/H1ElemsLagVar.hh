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
#include "FeBasis/FeNodal.hh"

namespace CoupledField {

  class FeH1LagrangeVar : public FeH1, public FeNodal {

  public:
    //! Constructor
    FeH1LagrangeVar();

    //! Copy Constructor
    FeH1LagrangeVar(const FeH1LagrangeVar& other);

      //! Destructor
    virtual ~FeH1LagrangeVar();

    //! Get number of shape functions for a given type (NODE/EDGE/FACE/ELEM)
    void GetNumFncs( StdVector<UInt>& numFcns,
                     EntityType fctEntityType,
                     UInt dof = 1 );

    //! Flag, if element has true nodal permutation
    virtual bool NeedsNodalPermutation() {
      return true;
    }
    
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

    //! Set the polynomial order of the element
    virtual void SetIsoOrder( UInt order ) = 0;
    
    //! \copydoc BaseFE::GetIsoOrder
    virtual UInt GetIsoOrder() const {
      return order_;
    }
    
    //! \copydoc BaseFE::GetMaxOrder
    virtual UInt GetMaxOrder() const {
      return order_;
    }

    //! @copydoc BaseFE::ComputeMonomialCoefficients
    //! Overloaded method for lagrange Elements
    virtual void ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C);

    //! Compare two element for equality (= same shape and approximation);
    bool operator==( const FeH1LagrangeVar& comp) const;

  protected:

    //! @copydoc FeNodal::SetFunctionsAtIp
    void SetFunctionsAtIp( const StdVector<LocPoint>& iPoints );

    //! @copydoc FeNodal::SetFunctionsAtIp
    void SetFunctionsAtIp( const std::map<Integer, LocPoint >& 
                           iPoints);

    //! @copydoc FeH1::CalcShFnc
    virtual void CalcShFnc( Vector<Double>& shape,
                            const Vector<Double>& point,
                            const Elem* ptElem,
                            UInt comp = 1 ) = 0;

    //! @copydoc FeH1::CalcLocDerivShFnc
    virtual void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                    const Vector<Double>& point,
                                    const Elem* ptElem,
                                    UInt comp = 1 ) = 0;
    


    //! Polynomial order of the finite element
    UInt order_;
  };

  //! Lagrangian line element of variable order
  class FeH1LagrangeLineVar : public FeH1LagrangeVar {

  public:

    //! Constructor
    FeH1LagrangeLineVar();

    //! Copy Constructor
    FeH1LagrangeLineVar(const FeH1LagrangeLineVar& other)
      :FeH1LagrangeVar(other) {
    }

    //! Destructor
    virtual ~FeH1LagrangeLineVar();

    virtual FeH1LagrangeLineVar* Clone(){
      return new FeH1LagrangeLineVar(*this);
    }

    //! Set the isotropic order of the Element. This methods gets overwritten 
    //! by the child classes to calculate the number of functions according to
    //! the given order
    //! \param order (input) The desired order of the element
    virtual void SetIsoOrder(UInt order);

  protected:

    //! @copydoc FeH1::CalcShFnc
    void CalcShFnc( Vector<Double>& shape,
                    const Vector<Double>& point,
                    const Elem* ptElem,
                    UInt comp = 1 );

    //! @copydoc FeH1::CalcLocDerivShFnc
    void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                            const Vector<Double>& point,
                            const Elem* ptElem,
                            UInt comp = 1 );

    //! Returns all coordinates of the given element in the correct ordering
    virtual void GetLocalDOFCoordinates(Matrix<Double> & coordMat);
  };
  
  //! Lagrangian quadrilateral element of variable order 
  class FeH1LagrangeQuadVar : public FeH1LagrangeVar {

  public:

    //! Constructor
    FeH1LagrangeQuadVar();

    //! Copy Constructor
    FeH1LagrangeQuadVar(const FeH1LagrangeQuadVar& other)
      :FeH1LagrangeVar(other) {
    }

    //! Destructor
    virtual ~FeH1LagrangeQuadVar();

    virtual FeH1LagrangeQuadVar* Clone(){
      return new FeH1LagrangeQuadVar(*this);
    }

    //! Set the isotropic order of the Element. This methods gets overwritten 
    //! by the child classes to calculate the number of functions according to
    //! the given order
    //! \param order (input) The desired order of the element
    virtual void SetIsoOrder(UInt order);

  protected:

    //! @copydoc FeH1::CalcShFnc
    void CalcShFnc( Vector<Double>& shape,
                    const Vector<Double>& point,
                    const Elem* ptElem,
                    UInt comp = 1 );

    //! @copydoc FeH1::CalcLocDerivShFnc
    void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                            const Vector<Double>& point,
                            const Elem* ptElem,
                            UInt comp = 1 );
    virtual void GetLocalDOFCoordinates(Matrix<Double> & coordMat);

    //! @copydoc FeH1::GetNumICModes
     virtual UInt GetNumICModes() {
    	 if ( order_ == 1 )
    		 return 2;
    	 else
    		 return 0;
     }

     //! @copydoc FeH1::CalcShFncICModes
     void CalcShFncICModes( Vector<Double>& shape,
                            const Vector<Double>& point,
                            const Elem* ptElem,
                            UInt comp = 1 );

     //! @copydoc FeH1::CalcLocDerivShFncICModes
     void CalcLocDerivShFncICModes( Matrix<Double> & deriv,
                                    const Vector<Double>& point,
                                    const Elem* ptElem,
                                    UInt comp = 1 );
  };
  
  //! Lagrangian hexahedral element of varaiable order
  class FeH1LagrangeHexVar : public FeH1LagrangeVar {

  public:

    //! Constructor
    FeH1LagrangeHexVar();

    //! Copy Constructor
    FeH1LagrangeHexVar(const FeH1LagrangeHexVar& other)
      :FeH1LagrangeVar(other) {
    }

    //! Destructor
    virtual ~FeH1LagrangeHexVar();

    virtual FeH1LagrangeHexVar* Clone(){
      return new FeH1LagrangeHexVar(*this);
    }

    //! Set the isotropic order of the Element. This methods gets overwritten 
    //! by the child classes to calculate the number of functions according to
    //! the given order
    //! \param order (input) The desired order of the element
    virtual void SetIsoOrder(UInt order);

  protected:

    //! @copydoc FeH1::CalcShFnc
    void CalcShFnc( Vector<Double>& shape,
                    const Vector<Double>& point,
                    const Elem* ptElem,
                    UInt comp = 1 );

    //! @copydoc FeH1::CalcLocDerivShFnc
    void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                            const Vector<Double>& point,
                            const Elem* ptElem,
                            UInt comp = 1 );

    virtual void GetLocalDOFCoordinates(Matrix<Double> & coordMat);
  };  
}
#endif  
