#ifndef FILE_SMOOTHINT
#define FILE_SMOOTHINT

#include <Elements/basefe.hh>
#include <Forms/bdbInt.hh>
#include <DataInOut/MaterialData.hh>

namespace CoupledField
{
  
  
  /// base class for calculation of linear elasticity
class SmoothInt : public BDBInt
{
public:
  /// Constructor
  SmoothInt(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  SmoothInt(MaterialData & matData);
  
  /// Destructor
virtual ~SmoothInt();
  
protected:    
  
  /// returns B - matrix for BDB
  virtual void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord);

  /// orientation of calculation plane in 2D 
  //  (especially important for anisotropic simulations)
  enum orientation2D {xy, xz, yz};

  orientation2D actOrientation;
};
  
  




  /// class for calculation of smoothanical plain strain state
class smoothPlainStrainInt : public SmoothInt
{  
public:
  /// Constructor
  smoothPlainStrainInt(BaseFE * aptelem, MaterialData & matDat);

  /// Constructor
  smoothPlainStrainInt(MaterialData & matDat);
  
  
  /// Deconstructor
  virtual ~smoothPlainStrainInt();
  
  
protected:
  
 /// calculate the data-matrix for 2D plain-strain
  virtual void calcDMat(Matrix<Double> & dMat, Integer ip, Matrix<Double> & ptCoord);

  /// returns dimension of D matrix
  virtual Integer getDimD(){return 3;};
  
  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs(){return 2;};
};


  /// class for calculation of smoothanical plain strain state
class smooth3DInt : public SmoothInt
{  
public:
  /// Constructor
  smooth3DInt(BaseFE * aptelem, MaterialData & matDat);

  /// Constructor
  smooth3DInt(MaterialData & matDat);
  
  /// Deconstructor
  virtual ~smooth3DInt();

  
protected:
  
  /// returns D - matrix for BDB
  virtual void calcDMat(Matrix<Double> & dMat, Integer ip, Matrix<Double> & ptCoord);

  /// returns dimension of D matrix
  virtual Integer getDimD(){return 6;};
  
  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs(){return 3;};  
};


} //end namespace

#endif // FILE_SMOOTHINT
