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
  FluidMechInt(Double density, Double kinematicViscosity, bool movingMesh, std::string stabilType="none");

  ///
  virtual ~FluidMechInt();

  /// Calculation of stiffmess matrix
  virtual void CalcElementMatrix( Matrix<Double>& stiffMat,EntityIterator& ent1,EntityIterator& ent2 ) {
    EXCEPTION( "FluidMechInt::CalcElementMatrix() not correctly overwritten!" );
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
  {EXCEPTION ("ResortElementVector is now in LinearForm");};

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

  inline void DyadicMult(Matrix<Double>& matrixRes, const Vector<Double>& vec1, \
      const Vector<Double>& vec2)
  {
    const UInt sizeRow = vec1.GetSize();
    const UInt sizeCol = vec2.GetSize();
    DyadicMult(matrixRes, vec1, vec2, sizeRow, sizeCol);
  }

  inline void DyadicMult(Matrix<Double>& matrixRes, const Vector<Double>& vec1, \
      const Vector<Double>& vec2, const UInt sizeRow, const UInt sizeCol)
  {
    Double* mtx_CurRow = matrixRes[0];
    for(UInt actRow=0; actRow < sizeRow; actRow++)
    {
      const Double& vec1_CurElem = vec1[actRow];
      for(UInt actCol=0; actCol < sizeCol; actCol++)
      {
        *mtx_CurRow = vec1_CurElem * vec2[actCol];
        ++mtx_CurRow;
      }
    }
  }

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

  bool movingMesh_;
  bool multAddOpt_;
  Double dt_;

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

protected:
  UInt spaceDim_;
  Matrix<Double> A_1_, A_2_;
  Matrix<Double> A_11_, A_12_;
  Matrix<Double> A_21_, A_22_;
  Matrix<Double> A_tmp_;
  Matrix<Double> xiDxDy_, xiDxxDyyDxy_;
  Vector<Double> DxXi_, DyXi_, xi_;
  Vector<Double> DxxXi_, DyyXi_;
  Vector<Double> auxI_, auxJ_;
  bool allocMem;

  Matrix<Double> elemResult_;
  Matrix<Double> locElemMat_UV_;
  Matrix<Double> locElemMat_PV_;
  Matrix<Double> locElemMat_UQ_;

  inline void AddSubMatrix(const Matrix<Double>& mat, \
      const Matrix<Double>& subMat, const UInt& startRow, const UInt& startCol)
  {
    const UInt& numRows = subMat.rows();
    const UInt& numCols = subMat.cols();
#ifdef CHECK_INITIALIZED
    if (subMat.rows() == 0 || subMat.cols() == 0 || mat.cols() == 0 || mat.rows() == 0 ) 
      EXCEPTION("undefined matrix" );
#endif
  
#ifdef CHECK_INDEX
    if ((subMat.rows() + startRow > mat.rows()) || (subMat.cols() + startCol > mat.cols()) )
      EXCEPTION("Submatrix to be read is to large!");
#endif

    for( UInt actRow=0; actRow < numRows; actRow++)
    {
      const UInt rowTmp = actRow + startRow;
      for( UInt actCol=0; actCol < numCols; actCol++)
      {
        mat[rowTmp][actCol + startCol] += subMat[actRow][actCol];
      }
    }
  }

  /* TODO: this should be done in initialisation phase, but the data is missing there */
  inline void doMemAlloc(const UInt& size)
  {
    allocMem = true;
    A_1_.Resize(size);
    A_2_.Resize(size);
    A_11_.Resize(size);
    A_12_.Resize(size);
    A_21_.Resize(size);
    A_22_.Resize(size);
    A_tmp_.Resize(size);
    DxXi_.Resize(size);
    DyXi_.Resize(size);
    DxxXi_.Resize(size);
    DyyXi_.Resize(size);
    locElemMat_UV_.Resize(spaceDim_ * size);
    locElemMat_PV_.Resize(spaceDim_ * size, size);
    locElemMat_UQ_.Resize(size, spaceDim_ * size);
  }
private:
};


}

#endif // FILE_LAPLACEINT
