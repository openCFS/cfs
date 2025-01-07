// =====================================================================================
// 
//       Filename:  curlOp.hh
// 
//    Description:  This operator implements the curl of the shape function
//                  Its operator-matrix has the structure
//                      / d N_1z/dy - d N_1y/dz ...\
//                  B = | d N_1x/dz - d N_1z/dx ...|
//                      \ d N_1y/dx - d N_1x/dy .../
//                      
//        Version:  1.0
//        Created:  10/28/2011 05:35:08 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#include "BaseBOperator.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "FeBasis/HCurl/HCurlElems.hh"

namespace CoupledField{

  template<class FE, UInt D, class TYPE>
  class CurlOperator : public BaseBOperator{
    public:
    
    CurlOperator();

    CurlOperator(const CurlOperator & other)
     : BaseBOperator(other){
    }

    virtual CurlOperator * Clone(){
      return new CurlOperator(*this);
    }

      ~CurlOperator();

      virtual void CalcOpMat(Matrix<Double> & bMat,
                             const LocPointMapped& lp, BaseFE* ptFe ){
        EXCEPTION("Curl Operator not implemented for general functions yet")
      }

      virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                       const LocPointMapped& lp, BaseFE* ptFe ){
        EXCEPTION("Curl Operator not implemented for general functions yet")
      }

      //avoid reimplementation of complex operator by making the bas class function
      //available, template caused
      using BaseBOperator::CalcOpMat;

      using BaseBOperator::CalcOpMatTransposed;

    protected:

  };

  template<class FE, UInt D, class TYPE>
   class CurlOperatorMag : public CurlOperator<FE,D,TYPE> {

   public:

    CurlOperatorMag(){
       this->name_ = "CurlOperatorMag";
     }

    CurlOperatorMag(const CurlOperatorMag & other)
      : CurlOperator<FE,D,TYPE>(other){
     }

     virtual CurlOperatorMag * Clone(){
       return new CurlOperatorMag(*this);
     }

     virtual ~CurlOperatorMag(){

     }

     virtual void CalcOpMat(Matrix<Double> & bMat,
                            const LocPointMapped& lp, BaseFE* ptFe );

     protected:

   };

  // Performing the calculation of (v x nabla x A)
  template<class FE, UInt D, class TYPE>
  void CurlOperatorMag<FE,D,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                         const LocPointMapped& lp,
                                         BaseFE* ptFe ){

    // bMatInitial containing the formfunctions in 2D(2x4)
    Matrix<Double> bMatInitial;
    CurlOperator<FE,D,TYPE>::CalcOpMat(bMatInitial,lp,ptFe);

    UInt numFncs = ptFe->GetNumFncs();

    // Initialize correct bMat
    if (bMatInitial.GetNumRows() == 2)
    {
      bMat.Resize(1, numFncs *1);
    }
    else
    {
      bMat.Resize(3, numFncs);
    }

    bMat.Init();

    //obtain external field velocity
    Vector<Double> myVec;
    int dof3 = 3;

    myVec.Resize(dof3);
    myVec.Init();
    this->coef_->GetVector(myVec,lp);
    Vector<Double> filler;
    filler.Resize(dof3);
    Vector<Double> buffer;
    buffer.Resize(dof3);

    // Calculate the crossproduct (v x nabla x Na). Vector filler is necessary because Na is a matrix of
    // all nodes (in 2D 2x4, in 3D 3x12), looping over the matrix brings the vector for the node or the edge,
    // buffer is the result of the crossproduct.

    for(unsigned int i=0; i< bMatInitial.GetNumCols(); ++i){
      filler.Init();
      buffer.Init();
//      std::cout << "erstes: myVec= " << myVec.ToString() << " filler= " << filler.ToString() << " buffer= " << buffer.ToString() << " bMatInitial= " << bMatInitial.ToString() << " bMat= " << bMat.ToString() << std::endl;
      for(unsigned int j=0; j< bMatInitial.GetNumRows(); ++j){
        filler[j] = bMatInitial[j][i];
      }
      myVec.CrossProduct(filler,buffer);
//      std::cout << "zweites: myVec= " << myVec.ToString() << " filler= " << filler.ToString() << " buffer= " << buffer.ToString() << " bMatInitial= " << bMatInitial.ToString() << " bMat= " << bMat.ToString() << std::endl;
      if (bMatInitial.GetNumRows() == 2)
      {
        double vel = buffer[2]; // because the result from crossproduct in 2D will be always in the third line
        bMat[0][i] = vel;
      }
      else
      {
        for(unsigned int u=0; u < buffer.GetSize(); ++u)
        {
          bMat[u][i] = buffer[u];
        }
      }
      //std::cout << "drittes: myVec= " << myVec.ToString() << " filler= " << filler.ToString() << " buffer= " << buffer.ToString() << " bMatInitial= " << bMatInitial.ToString() << " bMat= " << bMat.ToString() << std::endl;
    }

  }

  
  // ===================================
  //  CURL-OPERATOR FOR H-CURL ELEMENTS
  // ===================================
  //!Specialized class for HCurl elements
  template<UInt D, class TYPE >
  class CurlOperator<FeHCurl,D,TYPE> : public BaseBOperator{
    public:
    
    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 1;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = D; 
    //@}

    CurlOperator(){
      this->name_ = "CurlOperator";

    }

    CurlOperator(const CurlOperator & other)
     : BaseBOperator(other){
    }

    virtual CurlOperator * Clone(){
      return new CurlOperator(*this);
    }

    ~CurlOperator(){

    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe ){
      FeHCurl *fe = (static_cast<FeHCurl*>(ptFe));
      fe->GetCurlShFnc(bMat, lp, lp.shapeMap->GetElem(), 1);
    }

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe ){
      FeHCurl *fe = (static_cast<FeHCurl*>(ptFe));
      Matrix<Double> xiDx;
      fe->GetCurlShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
      xiDx.Transpose(bMat);
    }

    //avoid reimplementation of complex operator by making the bas class function
    //available, template caused
    using BaseBOperator::CalcOpMat;

    using BaseBOperator::CalcOpMatTransposed;

    // ===============
    //  QUERY METHODS
    // ===============
    //@{ \name Query Methods
    //! \copydoc BaseBOperator::GetDiffOrder
    virtual UInt GetDiffOrder() const {
      return ORDER_DIFF;
    }

    //! \copydoc BaseBOperator::GetDimDof()
    virtual UInt GetDimDof() const {
      return DIM_DOF;
    }

    //! \copydoc BaseBOperator::GetDimSpace()
    virtual UInt GetDimSpace() const {
      return DIM_SPACE;
    }

    //! \copydoc BaseBOperator::GetDimElem()
    virtual UInt GetDimElem() const {
      return DIM_ELEM;
    }

    //! \copydoc BaseBOperator::GetDimDMat()
    virtual UInt GetDimDMat() const {
      return DIM_D_MAT;
    }
    //@}
    
    protected:

  };

  // ===========================================================
  //  CURL-OPERATOR FOR H-CURL ELEMENTS (Scaled by Edge length)
  // ===========================================================
  
  //!Specialized class for curl elements
  template<UInt D, class TYPE>
  class ScaledByEdgeCurlOperator : public CurlOperator<FeHCurl,D,TYPE>{
  public:

    ScaledByEdgeCurlOperator(){

    }

    ScaledByEdgeCurlOperator(const ScaledByEdgeCurlOperator & other)
     : CurlOperator<FeHCurl,D,TYPE>(other){
    }

    virtual ScaledByEdgeCurlOperator * Clone(){
      return new ScaledByEdgeCurlOperator(*this);
    }

    ~ScaledByEdgeCurlOperator(){

    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe ){
      CurlOperator<FeHCurl,D,TYPE>::CalcOpMat(bMat,lp,ptFe);
      Double minE,maxE;
      lp.shapeMap->GetMaxMinEdgeLength(maxE,minE);
      bMat /= maxE;
    }

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe ){
      CurlOperator<FeHCurl,D,TYPE>::CalcOpMatTransposed(bMat,lp,ptFe);
      Double minE,maxE;
      lp.shapeMap->GetMaxMinEdgeLength(maxE,minE);
      bMat /= maxE;
    }

    //avoid reimplementation of complex operator by making the bas class function
    //available, template caused
    using CurlOperator<FeHCurl,D,TYPE>::CalcOpMat;

    using CurlOperator<FeHCurl,D,TYPE>::CalcOpMatTransposed;

  protected:

  };
  

  // ===========================================================
  //  CURL-OPERATOR FOR H-CURL ELEMENTS (Scaled by a mapping function)
  // ===========================================================

  //!Specialized class for curl elements
  template<UInt D, class TYPE>
  class ScaledCurlOperator : public CurlOperator<FeHCurl,D,TYPE>{
  public:
//TODO change the methods
    ScaledCurlOperator(){
      this->name_ = "ScaledCurlOperator";
    }

    ScaledCurlOperator(const ScaledCurlOperator & other)
     : CurlOperator<FeHCurl,D,TYPE>(other){
    }

    virtual ScaledCurlOperator * Clone(){
      return new ScaledCurlOperator(*this);
    }

    ~ScaledCurlOperator(){

    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe ){
      assert(this->coef_ != NULL);

      //TODO could be more difficult
      CurlOperator<FeHCurl,D,TYPE>::CalcOpMat(bMat,lp,ptFe);
      Vector<Double> coefs;
      this->coef_->GetVector(coefs,lp);
      for(UInt i=0;i<bMat.GetNumCols();++i){
        for(UInt d = 0;d<bMat.GetNumRows();++d){
          bMat[d][i] *= coefs[d];
        }
      }
    }

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe ){
      assert(this->coef_ != NULL);

      //TODO could be more difficult
      CurlOperator<FeHCurl,D,TYPE>::CalcOpMatTransposed(bMat,lp,ptFe);
      Vector<Double> coefs;
      this->coef_->GetVector(coefs,lp);
      for(UInt i=0;i<bMat.GetNumRows();++i){
        for(UInt d = 0;d<bMat.GetNumCols();++d){
          bMat[i][d] *= coefs[d];
        }
      }
    }


  protected:

  };

  // =====================================
  //  3D - CURL-OPERATOR FOR H-1 ELEMENTS
  // =====================================
  
  template<class TYPE>
  class CurlOperator<FeH1, 3, TYPE> : public BaseBOperator{
  public:

    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 3;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = 3;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = 3;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 3; 
    //@}
    
    CurlOperator(){
      this->name_ = "CurlOperator3D";

    }

    CurlOperator(const CurlOperator & other)
     : BaseBOperator(other){
    }

    virtual CurlOperator * Clone(){
      return new CurlOperator(*this);
    }

    ~CurlOperator(){

    }
    
    void CalcOpMat(Matrix<Double> & bMat,
                   const LocPointMapped& lp, BaseFE* ptFe ){

      UInt numFncs = ptFe->GetNumFncs();
      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( 3, numFncs * 3 );
      bMat.Init();

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      FeH1 *fe = (static_cast<FeH1*>(ptFe));
      fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );

      // see Kaltenbacher, 2nd edition, p. 124
      UInt i = 0;
      UInt pos = 0;
      for( i = 0, pos = 0; i < numFncs; ++i, pos+=3 ) {
        bMat[0][pos+1] = -xiDx[i][2];
        bMat[0][pos+2] =  xiDx[i][1];
      }

      for( i = 0, pos = 0; i < numFncs; ++i, pos+=3 ) {
        bMat[1][pos+0] =  xiDx[i][2];
        bMat[1][pos+2] = -xiDx[i][0];
      }

      for( i = 0, pos = 0; i < numFncs; ++i, pos+=3 ) {
        bMat[2][pos+0] = -xiDx[i][1];
        bMat[2][pos+1] =  xiDx[i][0];
      }
    }

    void CalcOpMatTransposed(Matrix<Double> & bMat,
                             const LocPointMapped& lp, BaseFE* ptFe ){
      UInt numFncs = ptFe->GetNumFncs();
      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( numFncs * 3, 3  );
      bMat.Init();

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      FeH1 *fe = (static_cast<FeH1*>(ptFe));
      fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );

      // see Kaltenbacher, 2nd edition, p. 124
      UInt i = 0;
      UInt pos = 0;
      for( i = 0, pos = 0; i < numFncs; ++i, pos+=3 ) {
        bMat[pos+0][1] =  xiDx[i][2];
        bMat[pos+0][2] = -xiDx[i][1];
      }

      for( i = 0, pos = 0; i < numFncs; ++i, pos+=3 ) {
        bMat[pos+1][0] = -xiDx[i][2];
        bMat[pos+1][2] =  xiDx[i][0];
      }

      for( i = 0, pos = 0; i < numFncs; ++i, pos+=3 ) {
        bMat[pos+2][1] = -xiDx[i][0];
        bMat[pos+2][0] =  xiDx[i][1];
      }
    }

    //avoid reimplementation of complex operator by making the bas class function
    //available, template caused
    using BaseBOperator::CalcOpMat;

    using BaseBOperator::CalcOpMatTransposed;

    // ===============
    //  QUERY METHODS
    // ===============
    //@{ \name Query Methods
    //! \copydoc BaseBOperator::GetDiffOrder
    virtual UInt GetDiffOrder() const {
      return ORDER_DIFF;
    }

    //! \copydoc BaseBOperator::GetDimDof()
    virtual UInt GetDimDof() const {
      return DIM_DOF;
    }

    //! \copydoc BaseBOperator::GetDimSpace()
    virtual UInt GetDimSpace() const {
      return DIM_SPACE;
    }

    //! \copydoc BaseBOperator::GetDimElem()
    virtual UInt GetDimElem() const {
      return DIM_ELEM;
    }

    //! \copydoc BaseBOperator::GetDimDMat()
    virtual UInt GetDimDMat() const {
      return DIM_D_MAT;
    }
    //@}

  protected:

  };

  // =====================================
  //  2D - CURL-OPERATOR FOR H-1 ELEMENTS
  // =====================================
  
  //! Curl operator for 2D fields
  
  //! Note: We assume, that the field itself is a vector with only 
  //! a z-component. 
  template<class TYPE>
  class CurlOperator<FeH1, 2, TYPE> : public BaseBOperator{
  public:

    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 1;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = 2;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = 2;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 2; 
    //@}
    
    CurlOperator(){
      this->name_ = "CurlOperator2D";

    }

    CurlOperator(const CurlOperator & other)
     : BaseBOperator(other){
    }

    virtual CurlOperator * Clone(){
      return new CurlOperator(*this);
    }

    ~CurlOperator(){

    }

    void CalcOpMat(Matrix<Double> & bMat,
                   const LocPointMapped& lpm, 
                   BaseFE* ptFe ){
      const UInt numFncs = ptFe->GetNumFncs();
      
      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( 2, numFncs );

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      FeH1 *fe = (static_cast<FeH1*>(ptFe));
      if (this->isSurfOpt_)
        fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
      else
        fe->GetGlobDerivShFnc( xiDx, lpm, lpm.shapeMap->GetElem() , 1 );

      UInt iFunc = 0;
      for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
        bMat[0][iFunc] =  xiDx[iFunc][1];
      }
      for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
        bMat[1][iFunc] = -xiDx[iFunc][0];
      }
    }

    void CalcOpMatTransposed(Matrix<Double> & bMat,
                             const LocPointMapped& lpm, BaseFE* ptFe ){
      
      const UInt numFncs = ptFe->GetNumFncs();

      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( numFncs, 2 );

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      FeH1 *fe = (static_cast<FeH1*>(ptFe));
      if (this->isSurfOpt_)
        fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
      else
        fe->GetGlobDerivShFnc( xiDx, lpm, lpm.shapeMap->GetElem() , 1 );
      
      UInt iFunc = 0;
      for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
        bMat[iFunc][0] =  xiDx[iFunc][1];
        bMat[iFunc][1] = -xiDx[iFunc][0];
      }
    }

    //avoid reimplementation of complex operator by making the base class function
    //available, template caused
    using BaseBOperator::CalcOpMat;

    using BaseBOperator::CalcOpMatTransposed;

    // ===============
    //  QUERY METHODS
    // ===============
    //@{ \name Query Methods
    //! \copydoc BaseBOperator::GetDiffOrder
    virtual UInt GetDiffOrder() const {
      return ORDER_DIFF;
    }

    //! \copydoc BaseBOperator::GetDimDof()
    virtual UInt GetDimDof() const {
      return DIM_DOF;
    }

    //! \copydoc BaseBOperator::GetDimSpace()
    virtual UInt GetDimSpace() const {
      return DIM_SPACE;
    }

    //! \copydoc BaseBOperator::GetDimElem()
    virtual UInt GetDimElem() const {
      return DIM_ELEM;
    }

    //! \copydoc BaseBOperator::GetDimDMat()
    virtual UInt GetDimDMat() const {
      return DIM_D_MAT;
    }
    //@}
    
  protected:
  };
  
  
  // ====================================================
  //  2D - AXISYMMETRIC - CURL-OPERATOR FOR H-1 ELEMENTS
  // ====================================================
   template<class TYPE>
   class CurlOperatorAxi : public BaseBOperator{
   public:
    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 1;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = 2;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = 2;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 2; 
    //@}
    
    CurlOperatorAxi(){
      this->name_ = "CurlOperatorAxi";

    }

    CurlOperatorAxi(const CurlOperatorAxi & other)
     : BaseBOperator(other){
    }

    virtual CurlOperatorAxi * Clone(){
      return new CurlOperatorAxi(*this);
    }

    ~CurlOperatorAxi(){

    }

    void CalcOpMat(Matrix<Double> & bMat,
                   const LocPointMapped& lpm, 
                   BaseFE* ptFe ){
      const UInt numFncs = ptFe->GetNumFncs();

      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( 2, numFncs );

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      FeH1 *fe = (static_cast<FeH1*>(ptFe));
      fe->GetGlobDerivShFnc( xiDx, lpm, lpm.shapeMap->GetElem() , 1 );

      // Calculate phi-phi component
      Vector<Double> shape;
      fe->GetShFnc( shape, lpm.lp, lpm.shapeMap->GetElem() );
      Vector<Double> globPoint;
      lpm.shapeMap->Local2Global(globPoint, lpm.lp);
      const Double oneOverR = 1.0 /  globPoint[0];

      UInt iFunc = 0;
      for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
        bMat[0][iFunc] = -xiDx[iFunc][1];
      }
      
      for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
        bMat[1][iFunc] =  xiDx[iFunc][0] + shape[iFunc] * oneOverR;
      }
    }

    void CalcOpMatTransposed(Matrix<Double> & bMat,
                             const LocPointMapped& lpm, 
                             BaseFE* ptFe ){
      
      const UInt numFncs = ptFe->GetNumFncs();

      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( numFncs, 2 );

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      FeH1 *fe = (static_cast<FeH1*>(ptFe));
      fe->GetGlobDerivShFnc( xiDx, lpm, lpm.shapeMap->GetElem() , 1 );
      
      // Calculate phi-phi component
      Vector<Double> shape;
      fe->GetShFnc( shape, lpm.lp, lpm.shapeMap->GetElem() );
      Vector<Double> globPoint;
      lpm.shapeMap->Local2Global(globPoint, lpm.lp);
      const Double oneOverR = 1.0 /  globPoint[0];

      UInt iFunc = 0;
      for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
        bMat[iFunc][0] = -xiDx[iFunc][1];
        bMat[iFunc][1] =  xiDx[iFunc][0] + shape[iFunc] * oneOverR;
      }
    }

    //avoid reimplementation of complex operator by making the base class function
    //available, template caused
    using BaseBOperator::CalcOpMat;

    using BaseBOperator::CalcOpMatTransposed;

    // ===============
    //  QUERY METHODS
    // ===============
    //@{ \name Query Methods
    //! \copydoc BaseBOperator::GetDiffOrder
    virtual UInt GetDiffOrder() const {
      return ORDER_DIFF;
    }

    //! \copydoc BaseBOperator::GetDimDof()
    virtual UInt GetDimDof() const {
      return DIM_DOF;
    }

    //! \copydoc BaseBOperator::GetDimSpace()
    virtual UInt GetDimSpace() const {
      return DIM_SPACE;
    }

    //! \copydoc BaseBOperator::GetDimElem()
    virtual UInt GetDimElem() const {
      return DIM_ELEM;
    }

    //! \copydoc BaseBOperator::GetDimDMat()
    virtual UInt GetDimDMat() const {
      return DIM_D_MAT;
    }
    //@}
    


  protected:  
   };
}
