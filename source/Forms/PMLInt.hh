// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PMLINT
#define FILE_PMLINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness/mass matrix for PML formulation

class PMLInt : public BaseForm
{
public:

  //! Constructor
  PMLInt(std::string type, Double factor, std::string dampingTypePML, 
         Double damp, bool axi=false);

  /// 
  virtual ~PMLInt();

  //! Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );

  //! set min/max of x,y,z coordinates form where PML starts and ends
  void SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer);
  
  void SetActElemSol(Matrix<Double>& disp){};


private:

  //! Calculation of stiffmess matrix
  void CalcElementMatrixStiff(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //! Calculation of mass matrix
  void CalcElementMatrixMass(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //! calculates position and values
  void ComputeFactorPML(Vector<Double>& factorsPML, Vector<Double>& coordAtIP);

  //! calculates the damping factor
  Double ComputeDampingFactor(Vector<Double>& pos, Directions dir);

  //! multiplicative factor for forms
  Double formsFactor_;

 //! type of bilinear form
  std::string formsType_;

  //! type of PML damping
  std::string dampingTypePML_;

  //! damping factor 
  Double dampingFactor_;

  //!layer thickness
  Matrix<Double> layerThickness_;

  //! coordinates for inner box at which PML starts
  Double minX_, maxX_, minY_, maxY_, minZ_, maxZ_;
};


}

#endif // FILE_LAPLACEINT
