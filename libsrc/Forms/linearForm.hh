#ifndef FILE_LINEARFORM_2
#define FILE_LINEARFORM_2

#include "baseForm.hh"

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




}

#endif // FILE_LINEARFORM
