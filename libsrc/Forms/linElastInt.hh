#ifndef FILE_LINELASTINT
#define FILE_LINELASTINT

#include <Elements/basefe.hh>
#include <Forms/bdbInt.hh>
#include <Materials/baseMaterial.hh>
#include <General/environment.hh>

namespace CoupledField {
  
  
  //! base class for calculation of linear elasticity
  class linElastInt : public BDBInt {

  public:

    //! Constructor
    linElastInt( BaseMaterial* matData, SubTensorType type = FULL );

    //! Destructor
    virtual ~linElastInt();

    //! Function for calculation bdb matrix 
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );

    //! just for computation of B - matrix
    void calcBMatOnly( Matrix<Double> &bMat, UInt ip,
		       BaseFE* elem, Matrix<Double> &ptCoord );

  protected:    

     //! calculates the material data 
    void calcDMat( Matrix<Double> &dMat );

     //! calculates the material data 
    void calcDMat( Matrix<Complex> &dMat );

    //! returns B - matrix for BDB
    virtual void calcBMat( Matrix<Double> &bMat, UInt ip,
                           Matrix<Double> &ptCoord );

    //! set dimensions
    virtual void SetDimensions(SubTensorType type);

    //! returns dimension of D matrix
    virtual UInt getDimD() {
      return dimD_; 
    };
    
    //! returns nr. of degrees of freedom
    virtual UInt getNrDofs() {
      return nrDofs_;
    };
    
  private:

    //dimension of Dmatrix
    UInt dimD_;
    
    //! number of degrees 
    UInt nrDofs_;

    //! subtype of tensor
    SubTensorType subTensorType_;
    
  };
  



} //end namespace

#endif // FILE_LINELASTINT
