#ifndef FILE_FLUIDMECHINT
#define FILE_FLUIDMECHINT

#include "Forms/baseForm.hh"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>

namespace CoupledField
{

  /// Class for calculation element stiffness matrix of a fluid mechanical operator
class FluidMechInt : public BaseForm
{
public:
  /// Constructor
  FluidMechInt(Double density, Double kinematicViscosity);

  /// 
  virtual ~FluidMechInt();

  /// Calculation of stiffmess matrix
  virtual void CalcElementMatrix( Matrix<Double>& stiffMat,EntityIterator& ent1,EntityIterator& ent2 ) {
    Error( "FluidMechInt::CalcElementMatrix() not correctly overwritten!", __FILE__, __LINE__);
  };


  virtual void SetActElemSol(Matrix<Double>& elemSol){};

//   void SetGridSolution( NodeStoreSol<Double>& gridSol ) {
//     gridSol_ = &gridSol; }
  
protected: 

  void ResortElementMatrix(Matrix<Double> & sortedElemMat, 
                           const Matrix<Double> & elemMat,
                           const UInt & nrOfNodes, 
                           const UInt & dofPerNode);
  void ColResortElementMatrix(Matrix<Double> & sortedElemMat, 
                              const Matrix<Double> & elemMat,
                              const UInt & nrOfFctRow, 
                              const UInt & nrOfFctCol, 
                              const UInt & dofPerNodeRow,
                              const UInt & dofPerNodeCol);
  void RowResortElementMatrix(Matrix<Double> & sortedElemMat, 
                              const Matrix<Double> & elemMat,
                              const UInt & nrOfFctRow, 
                              const UInt & nrOfFctCol, 
                              const UInt & dofPerNodeRow,
                              const UInt & dofPerNodeCol);
  
  void ResortElementVector(Vector<Double> & sortedElemVec, 
                           const Vector<Double> & elemVec,
                           const UInt & nrOfNodes, 
                           const UInt & dofPerNode)
  {Error("ResortElementVector is now in LinearForm",__FILE__,__LINE__);};
  
  void ComputeTaus(Double & tau_mp,
                   Double & tau_mu,
                   Double & tau_c,
                   const Double & VL2, 
                   const Double & VL2AtIp, 
                   const Double & VMax, 
                   const Double & hk, 
                   const Double & kinematicViscosity, 
                   const Double & lambda_k,
                   const UInt & nrFncs);

  //  void PrepareFluidMechIntegrator(const UInt nrFncs, const UInt nrIntPts);

  //Double ComputeTauM();
  //Double ComputeTauP();
  

  /// multiplicative value for laplace integration 
  Double density_;
  Double kinematicViscosity_;

  Double convStabSign_, viscStabSign_, presStabSign_, contiStabSign_, penaltyStabSign_;
  Double presTermSign_, contiTermSign_, penaltyTermSign_;
  Double penalty_;

//  Double tau_mpFac_, tau_muFac_, tau_cFac_;
  std::string stabilType_;
  Double velNormMin_, velNormMax_;

//  Double CFLMin_, CFLMax_;
//  Double PeeMin_, PeeMax_;

  std::string stabParamEsti_; //, movingMeshStr_;

  bool eigenValue_;
  bool movingMesh_;
  bool multAddOpt_;
  Double dt_;

  bool resortedWay_;

  //! grid solution vector
  //NodeStoreSol<Double>* gridSol_;

  std::vector<std::vector <Integer> > lookuptable;
  //Vector<Vector <Integer> > lookuptable;

//   UInt elemNumber_;
//   Matrix<Double> elemResult_, gridElemResult_;
//   Vector<Double>  Vx_, Vy_, Vz_;  
//   Double  VxAtIP_, VyAtIP_, VzAtIP_, VL2_, VL2AtIP_, VMax_;
//   Double tau_m_, tau_mu_, tau_c_;
//   Double lambda_k_, A_elem_, h_k_;

//   bool computeTaus_;

//   Double jacDet_, multAux_;
//   UInt N_; // DOFs per Node z.B.:(Ux,Uy)

//   Matrix<Double> locElemMat_;

//   Matrix<Double> A1_;
//   Matrix<Double> A_a_; 

//   Matrix<Double> xiDxDy_;
//   Matrix<Double> xiDxxDyyDxy_;

//   Vector<Double>  xiDxx_, xiDyy_;//, xiDxy;
//   Vector<Double>  xiDx_, xiDy_;
//   Vector<Double>  xi_;
//   Vector<Double>  auxI_, auxJ_;

private:
};


}

#endif // FILE_LAPLACEINT
