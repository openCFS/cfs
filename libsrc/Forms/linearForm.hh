#ifndef FILE_LINEARFORM_2
#define FILE_LINEARFORM_2

#include "baseForm.hh"
#include "Forms/nLinElastInt.hh"
#include "Utils/ApproxData.hh"

namespace CoupledField
{

/// base class class for calculation right hand side
class LinearForm : public BaseForm
{
public:
  ///
  LinearForm(BaseFE * aptelem);

  ///
  LinearForm();

  ///
  virtual ~LinearForm();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & result);


};



// =============================================================================
// edge integration
// =============================================================================


/// class for calculation of right hand side of edge elements
class LinearEdgeInt : public LinearForm
{
public:
  ///
  LinearEdgeInt(BaseFE * aptelem, Double val, Integer direction, 
		Vector<Double> * coilMidPoint = NULL);

  ///
  virtual ~LinearEdgeInt();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & result);

private:
  /// source factor
  Double val_;
  
  /// direction of source
  /*! 1: x-direction, 2: y-direction, 3: z-direction
   */
  Integer direction_;

  /// midpoint of coil (needed for circular coils to calculate the current dirction)
  Vector<Double> * coilMidPt_;  
};




// =============================================================================
// volume source integration
// =============================================================================


/// class for calculation of right hand side of an volume source
class VolumeSrcInt : public LinearForm
{
public:
  ///
  VolumeSrcInt(Double val, Boolean isaxi);

  ///
  virtual ~VolumeSrcInt();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & result);

private:
  /// source factor
  Double val_;
};


// =============================================================================
// permanent magnets in 2D
// =============================================================================


/// class for calculation of right hand side of a permanent magnet
class MagPerm2DInt : public LinearForm
{
public:
  ///
  MagPerm2DInt(Vector<Double> val, Double rel, Boolean isaxi);

  ///
  virtual ~MagPerm2DInt();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & result);

private:
  //! magnetization
  Vector<Double> perm_;

  //!reluctivity
  Double reluctivity_;
};


// =============================================================================
// nonlinear magnetics
// =============================================================================

/// class for calculation of right-hand-side of nonlinear magnetics
class nLinMagNode2D_linFormInt : public LinearForm
{
public:
  /// constructor
   nLinMagNode2D_linFormInt(BaseFE * aptelem, MaterialData & matData, Boolean axi=FALSE);

  /// constructor
   nLinMagNode2D_linFormInt(ApproxData *nlinFnc, Double startVal, Boolean axi=FALSE);

  /// constructor for linear subdomain
   nLinMagNode2D_linFormInt(Double startVal, Boolean axi=FALSE);

  /// destructor
  virtual ~ nLinMagNode2D_linFormInt();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & result);

 /// in nonlinear calculations, the actual magnetic vector potential of the element is needed
  virtual void SetActElemSol(Matrix<Double>& magPot) 
  { magPotinMatrix_ = magPot; magPot.ConvertToVec_AppendRows(magPot_);};
  
protected:
  /// material data
  MaterialData matData_;

  Double startmatVal_;
  ApproxData *nlinFnc_;
  Vector<Double> magPot_;
  Matrix<Double> magPotinMatrix_;
};


// =============================================================================
// nonlinear mechanics
// =============================================================================


/// class for calculation of right-hand-side of nonlinear mechanics
class nLinMech_linFormInt : public LinearForm
{
public:
  /// constructor
  nLinMech_linFormInt(BaseFE * aptelem, MaterialData & matData, Boolean axi=FALSE);

  /// constructor
  nLinMech_linFormInt(MaterialData & matData, Boolean axi=FALSE);

  /// destructor
  virtual ~nLinMech_linFormInt();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & result);

 /// in nonlinear calculations, the actual displacement of the element is needed
  /*!
  \param disp (input) Matrix with displacement d of all nodes of actual element
  \f[ \left( \begin{array}{ccc} 
             d_{x1} &  d_{x2} &  d_{x3} \\
             d_{y1} &  d_{y2} &  d_{y3} \\
	     \end{array}\right) \f]	    
  */
  void setActElemDispl(Matrix<Double>& disp) {elemDisp_ = disp;};  


  /// in nonlinear calculations, the actual displacement of the element is needed
  /*!
  \param disp (input) Matrix with displacement d of all nodes of actual element
  \f[ \left( \begin{array}{ccc} 
             d_{x1} &  d_{x2} &  d_{x3} \\
             d_{y1} &  d_{y2} &  d_{y3} \\
	     \end{array}\right) \f]	    
  */
  virtual void SetActElemSol(Matrix<Double>& disp) {elemDisp_ = disp;};
 
  
protected:
  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs(){return 3;};

  /// material data
  MaterialData matData_;

  /// displacement of all nodes of actual element
  Matrix<Double> elemDisp_;
};





// =============================================================================
// prestress
// =============================================================================



/// class for calculation of right-hand-side of prestress
class PreStressLinFormInt : public nLinMech_linFormInt
{
public:
  /// constructor
  PreStressLinFormInt(BaseFE * aptelem, 
		      MaterialData & matData, 
		      Double aPreStressVal, 
		      Directions stressDir);
  

  /// destructor
  virtual ~PreStressLinFormInt();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, 
			      Vector<Double> & result);

private: 
  ///
  Double preStressVal_;

  ///
  Directions preStressDir_;
};


// =============================================================================
// pressure load
// =============================================================================


/// class for surface integration
class PressureLinForm : public LinearForm
{
public:
  ///
  PressureLinForm(Double aVal, Boolean isaxi);

  ///
  virtual ~PressureLinForm();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, 
			      Vector<Double> & elemVec);

  virtual void SetMultiplier(Double mult){multiplier_ = mult;};
  
protected:

private:
  /// factor of load
  Double multiplier_;

};



/// class for calculation of right hand side for flownoise problem
class LinearFlowNoiseInt : public LinearForm
{
public:
  ///
  LinearFlowNoiseInt(BaseFE * aptelem);

  ///
  virtual ~LinearFlowNoiseInt();

  /// Calculation of vector of right hand side for the surface elements on the obstacle (dipole)
  void CalcElemVector4Dip(Matrix<Double>& ptCoord, 
			  const StdVector<Integer> & connecth, 
			  Vector<Double> & Result, 
			  const Vector<Double> gradN_x_P);

  /// Calculation of vector of right hand side given from quadrupole contribution
  void CalcElemVector4Quad(Matrix<Double>& ptCoord, 
			   const StdVector<Integer> & connecth,
			   const Matrix<Double> & FlowData, 
			   Vector<Double> & Result);

  /// Extraction of element velocity values from total flowdata matrix to a matrix (connecth, dim)
  void GetQttiesOfElement(Matrix<Double>& elVec, 
			  const Matrix<Double>& FlowData,
			  const StdVector<Integer>& connecth, 
			  Integer matrixRow);
  
  
private:

};


/// class for calculation of right hand side for recovery procedure
class RHSForRecoveryProcedure : public LinearForm
{
public:
  ///
  RHSForRecoveryProcedure(BaseFE * aptelem);

  ///
  virtual ~RHSForRecoveryProcedure();

  // Calculation of vector of right hand side
  /*
    /\f[
    \int \phi \frac{\partial u^{FEM}}{\partial x_i}
    \f]
    \param ptCoord (input) Matrix with coordinates of the element
    \param fncNodesElem (input) value of solution at nodes of element
    \param aComponent (input) number of the gradient's component 
    \param elemVec (output) vector with result
  */
  void CalcElemVectorRHSForSPR(Matrix<Double>& ptCoord,
			  Vector<Double> & fncNodesElem,
			  const Integer aComponent,
			  Vector<Double> & elemVec);

private:
 
};

}

#endif // FILE_LINEARFORM
