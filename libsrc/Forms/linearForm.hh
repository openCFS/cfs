#ifndef FILE_LINEARFORM_2
#define FILE_LINEARFORM_2

#include "baseForm.hh"
#include <Forms/nLinElastInt.hh>

namespace CoupledField
{

/// base class class for calculation right hand side
class LinearForm : public BaseForm
{
public:
  ///
  LinearForm(BaseFE * aptelem);

  ///
  virtual ~LinearForm();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & result);
  
};



/// class for calculation of right hand side of edge elements
class LinearEdgeInt : public LinearForm
{
public:
  ///
  LinearEdgeInt(BaseFE * aptelem, Double val, Integer direction, 
		std::vector<Double> * coilMidPoint = NULL);

  ///
  virtual ~LinearEdgeInt();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & result);

private:
  /// source factor
  Double val_;
  
  /// direction of source
  /*! 1: x-direction, 2: y-direction, 3: z-direction
   */
  Integer direction_;

  /// midpoint of coil (needed for circular coils to calculate the current dirction)
  std::vector<Double> * coilMidPt_;  
};




/// class for calculation of right-hand-side of nonlinear mechanics
class nLinMech_linFormInt : public LinearForm
{
public:
  /// constructor
  nLinMech_linFormInt(BaseFE * aptelem, MaterialData & matData);

  /// destructor
  virtual ~nLinMech_linFormInt();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & result);

 /// in nonlinear calculations, the actual displacement of the element is needed
  /*!
  \param disp (input) Matrix with displacement d of all nodes of actual element
  \f[ \left( \begin{array}{ccc} 
             d_{x1} &  d_{x2} &  d_{x3} \\
             d_{y1} &  d_{y2} &  d_{y3} \\
	     \end{array}\right) \f]	    
  */
  void setActElemDispl(Matrix<Double>& disp) {elemDisp_ = disp;};  


private:
  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs(){return 3;};

  /// material data
  MaterialData matData_;

  /// displacement of all nodes of actual element
  Matrix<Double> elemDisp_;
};




}

#endif // FILE_LINEARFORM
