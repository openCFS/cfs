#ifndef FILE_BASEFORM_
#define FILE_BASEFORM_

#include "Utils/StdVector.hh"
#include "Utils/vector.hh"
#include "Elements/basefe.hh"
#include "DataInOut/MaterialData.hh"

namespace CoupledField
{

  //! Base class for calculation of bdb element matrices
  class BaseForm 
  {
  public:

    //! Constructor
    BaseForm(BaseFE * aptelem, MaterialData & matData);

    //! Constructor
    BaseForm(MaterialData & matData);

    //! Constructor
    BaseForm(BaseFE * aptelem);

    //! Constructor
    BaseForm();

    //! Deconstructor
    virtual ~BaseForm();

    //! Virtual function, implemented in derived classes
    virtual void CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double> & StiffMat);

    //! Virtual function, implemented in derived classes
    virtual void CalcComplexElementMatrix(Matrix<Double>& ptCoord,
                                          Matrix<Complex> & StiffMat,
                                          Double & beta, Double & omega);
  
    /// Calculation of vector of right hand side 
    virtual void CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & result)
    {Error("CalcElemVector not implemented!",__FILE__,__LINE__);};

    virtual void CalcElemVector4Dip(Matrix<Double>& ptCoord, 
                                    const StdVector<Integer> & connecth, 
                                    Vector<Double> & Result, 
                                    const Vector<Double> gradN_x_P)
    { Error(" CalcElemVector4Dip is not implemented for this class",__FILE__,__LINE__);};

    /// Calculation of vector of right hand side given from quadrupole contribution
    virtual void CalcElemVector4Quad(Matrix<Double>& ptCoord,
                                     const StdVector<Integer> & connecth,
                                     const Matrix<Double> & FlowData, 
                                     Vector<Double> & Result)
    { Error(" CalcElemVector4Quad is not implemented for this class",__FILE__,__LINE__);};

    /// Extraction of element velocity values from total flowdata matrix to a matrix (connecth, dim)
    virtual void GetQttiesOfElement(Matrix<Double>& elVec,
                                    const Matrix<Double>& FlowData,
                                    const StdVector<Integer>& connecth, 
                                    Integer matrixRow)
    { Error(" GetQttiesOfElement is not implemented for this class",__FILE__,__LINE__);};

    //! Prints the bilinear form
    virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

    //
    virtual void SetNonLinMethod(std::string atype) {;};

    //! needed for fractional damping model
    virtual void SetFracDamping() {;};

    //! needed for fractional damping model
    virtual void UnsetFracDamping() {;};

    //! needed for fractional damping model
    virtual Boolean IsFracDamping()
    {return isFracDamping_;};

    //! needed for fractional damping model
    virtual void SetRaylDamping() 
    { isRaylDamping_ = TRUE; };


    //! needed for fractional damping model
    virtual Boolean IsRaylDamping()
    {return isRaylDamping_;};

    //! set additional multiplicative factor for matrix
    virtual void SetFactor(Double factor) {;};

    //! set second multiplicative factor for matrix
    virtual void SetSecondFactor(Double factor) {;};

    //! sets pointer to actual element
    void SetElemPtr(BaseFE * elemPtr){ptelem = elemPtr;};

    //! sets pointer to actual material
    void SetMaterial(MaterialData * matPtr){ptMaterial = matPtr;};

    MaterialData * GetMaterial(){return ptMaterial;};

    //! get base type of biliearform: STIFFNESS, DAMPING, MASS
    virtual FEMatrixType GetBaseType() 
    { return baseType_; };

    //! sets actual element solution
    virtual void SetActElemSol(Matrix<Double>& disp)
    {Error("SetActElemSol not implemented!",__FILE__,__LINE__);};

    //! sets actual element solution
    virtual void SetActElemSol(CFSMatrix & disp)
    {Error("SetActElemSol not implemented!",__FILE__,__LINE__);};

    //! sets actual element solution
    virtual void SetActElemSol(Vector<Double>& disp)
    {Error("SetActElemSol not implemented!",__FILE__,__LINE__);};

    //! sets actual first time derivative of element solution
    virtual void SetActElemSolDeriv1(Matrix<Double>& disp)
    {Error("SetActElemSol not implemented!",__FILE__,__LINE__);};

    //! sets actual first time derivative of element solution
    virtual void SetActElemSolDeriv1(Vector<Double>& disp)
    {Error("SetActElemSolDeriv1 not implemented!",__FILE__,__LINE__);};

    //! sets actual 2nd time derivative of element solution
    virtual void SetActElemSolDeriv2(Vector<Double>& disp)
    {Error("SetActElemSolDeriv2 not implemented!",__FILE__,__LINE__);};

    //! reads the values y(x) out of the file with name fncName  
    void ReadNlinFunc(std::string fncName, Vector<double> &xval, Vector<Double> &yval)
    {Error("ReadNlinFunc not implemented!");};

    //! set the integration point
    void SetIntPoint(Vector<Double> point)
    { intPoint_ = point; isSetIntPoint_ = TRUE;};

    //!
    void SetDofZero(Integer posdof)
    {dofzero_ = posdof; };

    //!
    void SetPiezoMaterialType(piezoMaterialType & pMatType)
    { piezoMatType_=pMatType; };

    //!
    void SetMaterialArray(Matrix<Double>* mat)
    { materialArray_ = mat; };

    //!
    void SetSubdomain(Integer sd)
    {actSD_ = sd; };

    //!
    void SetElemNr(Integer nr)
    {actElemNr_ = nr; };

  protected:


    BaseFE  * ptelem;   //!< pointer to base element
    MaterialData * ptMaterial ;   //!< pointer to material data
    Boolean isaxi_;  //!< true for axisymmetric setup

    //
    Vector<Double> intPoint_;
    //
    Boolean isSetIntPoint_;

    piezoMaterialType piezoMatType_;     //! default = realMaterialParamter, piezoMatType_ = imagMaterialParamter if we consider complex-valued material Paramter;

    Boolean isFracDamping_;   //!< if true Assemble::AssembleMatrices will retrieve an additional multiplicative factor
    Boolean isRaylDamping_;   //!< if true Assemble::AssembleMatrices will retrieve an additional multiplicative factor
    Integer dofzero_;   //!< for multidof-handling, where one dof is zero (e.g. piezoelectric PDE)

    FEMatrixType baseType_;  // base type: STIFFNESS, DAMPING, MASS

    Matrix<Double>* materialArray_;

    Integer actSD_;
    Integer actElemNr_;

  };

} //end namespace

#endif // FILE_BASEFORM
