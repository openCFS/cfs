#ifndef FILE_NLINELASTINT_03
#define FILE_NLINELASTINT_03

#include <Elements/basefe.hh>
#include <DataInOut/MaterialData.hh>

#include <Forms/linElastInt.hh>

namespace CoupledField
{
  

  /// base class for calculation of nonlinear linear elasticity
class nLinElastInt : public linElastInt
{
public:
  /// Constructor
  nLinElastInt(BaseFE * aptelem, MaterialData & matData);
  
  /// Destructor
  virtual ~nLinElastInt();  

  /// in nonlinear calculations, the actual displacement of the element is needed
  /*!
  \param disp (input) Matrix with displacement d of all nodes of actual element
  \f[ \left( \begin{array}{ccc} 
             d_{x1} &  d_{x2} &  d_{x3} \\
             d_{y1} &  d_{y2} &  d_{y3} \\
	     \end{array}\right) \f]	    
  */
  void setActElemDispl(Matrix<Double>& disp) {elemDisp_ = disp;};  


protected:    

  /// returns B - matrix for BDB
  virtual void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord);


  /// displacement of all nodes of actual element
  Matrix<Double> elemDisp_;
};
  
  
 


  /// class for calculation of 3d nonlinear linear elasticity
// derived from nLinElastInt with 3d material 
// matrix and basic nonlinearity methods
class nLinMech3dInt : public nLinElastInt
{
public:


  /// Constructor
  nLinMech3dInt(BaseFE * aptelem, MaterialData & matData);
  
  /// Destructor
  virtual ~nLinMech3dInt();  
  
protected:  
  /// returns D - matrix for BDB
  virtual void calcDMat(Matrix<Double> & dMat);

  /// returns dimension of D matrix
  virtual Integer getDimD(){return 6;};
  
  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs(){return 3;};  
};
  


  
} //end namespace

#endif // FILE_XXX


