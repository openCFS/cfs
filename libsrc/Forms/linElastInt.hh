#ifndef FILE_LINELASTINT
#define FILE_LINELASTINT

#include <Elements/basefe.hh>
#include <Forms/bdbInt.hh>
#include <DataInOut/MaterialData.hh>

namespace CoupledField
{
  
  
  /// base class for calculation of linear elasticity
class linElastInt : public BDBInt
{
public:
  /// Constructor
linElastInt(BaseFE * aptelem, MaterialData & matData);
  
  /// Destructor
virtual ~linElastInt();
  
protected:    
  
  /// returns B - matrix for BDB
  virtual void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord);

  /// orientation of calculation plane in 2D 
  //  (especially important for anisotropic simulations)
  enum orientation2D {xy, xz, yz};

  orientation2D actOrientation;
};
  
  




  /// class for calculation of mechanical plain strain state
class mechPlainStrainInt : public linElastInt
{  
public:
  /// Constructor
  mechPlainStrainInt(BaseFE * aptelem, MaterialData & matDat);
  
  
  /// Deconstructor
  virtual ~mechPlainStrainInt();
  
  
protected:
  
 /// calculate the data-matrix for 2D plain-strain
  virtual void calcDMat(Matrix<Double> & dMat);

  /// returns dimension of D matrix
  virtual Integer getDimD(){return 3;};
  
  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs(){return 2;};
};


  /// class for calculation of mechanical plain strain state
class mech3DInt : public linElastInt
{  
public:
  /// Constructor
  mech3DInt(BaseFE * aptelem, MaterialData & matDat);
  
  /// Deconstructor
  virtual ~mech3DInt();

  
protected:
  
  /// returns D - matrix for BDB
  virtual void calcDMat(Matrix<Double> & dMat);

  /// returns dimension of D matrix
  virtual Integer getDimD(){return 6;};
  
  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs(){return 3;};  
};


} //end namespace

#endif // FILE_BASERFORM
