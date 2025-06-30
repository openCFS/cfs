#ifndef FILE_COEFFUNCTION_SURF_HH
#define FILE_COEFFUNCTION_SURF_HH


#include "CoefFunction.hh"
#include "CoefFunctionMulti.hh"
namespace CoupledField {

// forward class declaration
struct ResultInfo;

//! Evaluates a coefficient function on a surface 

//! This class represents coefficient functions, which are defined just on a
//! surface. Thus, it internally holds pointers to coefficient functions for
//! each region. If queried for a scalar / vector / tensor entry of a given
//! surface element, it looks for the "correct" volume coefficient function
//! and returns it. In addition, it can also perform a normal-mapping of
//! the coefficient function.
class CoefFunctionSurf : public CoefFunction {
public:

  //! Constructor

  //! Constructor for the class
  //! \param mapNormal If true, only the normal component w.r.t. to the
  //!                  surface element is taken into account.By default,
  //!                  the normal direction points OUT of the related volumes.
  //! \param factor Additional scaling factor 
  //! \param surfInfo Result info object for surface result

  CoefFunctionSurf( bool mapNormal,
                    Double factor = 1.0,
                    shared_ptr<ResultInfo> surfInfo =  shared_ptr<ResultInfo>());

  //! Destructor
  virtual ~CoefFunctionSurf();

  virtual string GetName() const { return "CoefFunctionSurf"; }

  //! Set single volume coefficient function
  virtual void AddVolumeCoef( RegionIdType, PtrCoefFct );
  
  //! Pass volume coefficients
  virtual void SetVolumeCoefs( std::map<RegionIdType, PtrCoefFct> coefs );

  //! \copydoc CoefFunction::GetTensor
  virtual void GetTensor(Matrix<Double>& coefMat,
                 const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetTensor
  virtual void GetTensor(Matrix<Complex>& coefMat,
                   const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<Double>& coefVec,
                 const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<Complex>& coefVec,
                 const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetScalar
  virtual void GetScalar(Complex& coefScalar,
                   const LocPointMapped& lpm );

  // Workaround to use templated functions
  //! \copydoc CoefFunction::GetVector
  template<typename TYPE> void GetVector_(Vector<TYPE>& coefVec, const LocPointMapped& lpm );

    // Workaround to use templated functions
  //! \copydoc CoefFunction::GetVector
  template<typename TYPE> void GetScalar_(TYPE& coefScalar, const LocPointMapped& lpm );

  

  virtual void SetSurfInfo(LocPointMapped &surfLpm);

  //! \copydoc CoefFunction::GetVecSize
  UInt GetVecSize() const;

  //! Set the dimension of the physical vector returned by the coef function (needed in CoefExpressions)
  void SetDim(UInt dim){
    dim_ = dim;
  }

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const;

  //! \copydoc CoefFunction::ToString
  std::string ToString() const;

  // ===========================
  //  NORMAL MAPPING OPERATIONS
  // ===========================
  //@{ \name Functions for normal mapping
  
  //! Mapping operation for scalars
  template<typename TYPE>
  static void MapScalNormal( Vector<TYPE>& ret, TYPE& scal,
                            const Vector<Double>& normal );

  //! Mapping operation for vectors
  template<typename TYPE>
  static void MapVecNormal( TYPE& ret, const Vector<TYPE>& vec, 
                            const Vector<Double>& normal );
  
  //! Mapping operation for tensors in Voigt notation
  template<typename TYPE>
  static void MapTensorNormal( Vector<TYPE>& ret, const Vector<TYPE>& tensor,
                               const Vector<Double>& normal );
  
  //@}
  

  //! Map with CoefFunctions for each region
  std::map<RegionIdType, PtrCoefFct> coefs_;
  
  //! Set with all regionIdTypes
  std::set<RegionIdType> regions_;

  //! Flag, if normal mapping should be performed
  bool mapNormal_;
  
  //! Factor in case of surface mapping
  Double factor_;

  //! information
  SolutionType resultType_;

  //! dimension of space, needed to return sipe of physical vectors
  UInt dim_;
};

//! This class represents coefficient functions, which are defined just on a
//! surface and computes the force defined by Maxwell's stress tensor
//! It#s derived from CoefFunctionSurf
class CoefFunctionSurfMaxwell : public CoefFunctionSurf {
public:

  //! Constructor

  //! Constructor for the class
  //! \param mapNormal If true, only the normal component w.r.t. to the
  //!                  surface element is taken into account.By default,
  //!                  the normal direction points OUT of the related volumes.
  //! \param material parameter
  //! \param factor Additional scaling factor
  //! \param surfInfo Result info object for surface result

  CoefFunctionSurfMaxwell( bool mapNormal,
		                       std::map<SolutionType, shared_ptr<CoefFunctionMulti> > matCoefs,
                           Grid* ptrid, Double factor = 1.0,
                           shared_ptr<ResultInfo> surfInfo =  shared_ptr<ResultInfo>());

  //! Destructor
  virtual ~CoefFunctionSurfMaxwell();

  virtual string GetName() const { return "CoefFunctionSurfMaxwell"; }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat,
                 const LocPointMapped& lpm ) {
	  EXCEPTION("CoefFunctionSurfMaxwell:GetTensor not implemented");
  }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Complex>& coefMat,
		  const LocPointMapped& lpm ) {
  	  EXCEPTION("CoefFunctionSurfMaxwell:GetTensor not implemented");
  }

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<Double>& coefVec,
                 const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<Complex>& coefVec,
                 const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(Double& coefScalar,
                 const LocPointMapped& lpm ) {
  	  EXCEPTION("CoefFunctionSurfMaxwell:GetScalar not implemented");
  }

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(Complex& coefScalar,
		  const LocPointMapped& lpm ) {
  	  EXCEPTION("CoefFunctionSurfMaxwell:GetScalar not implemented");
  };


private:

  //! coef-function as defined in PDE
  std::map<SolutionType, shared_ptr<CoefFunctionMulti> > matCoef_;

  //! pointer to the grid
  Grid* ptGrid_;
};


//! This class represents coefficient functions, which are defined just on a
//! surface and computes the force defined by virtual work principle
//! It's derived from CoefFunctionSurf
class CoefFunctionSurfVWP : public CoefFunctionSurf {
public:

  //! Constructor

  //! Constructor for the class
  //! \param mapNormal If true, only the normal component w.r.t. to the
  //!                  surface element is taken into account.By default,
  //!                  the normal direction points OUT of the related volumes.
  //! \param material parameter
  //! \param factor Additional scaling factor
  //! \param surfInfo Result info object for surface result

  CoefFunctionSurfVWP( bool mapNormal,
	                   std::map<SolutionType, shared_ptr<CoefFunctionMulti> > matCoefs,
                       Double factor = 1.0,
					   shared_ptr<ResultInfo> surfInfo =  shared_ptr<ResultInfo>());

  //! Destructor
  virtual ~CoefFunctionSurfVWP();

  virtual string GetName() const { return "CoefFunctionSurfVWP"; }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat,
                 const LocPointMapped& lpm ) {
	  EXCEPTION("CoefFunctionSurfVWP:GetTensor not implemented");
  }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Complex>& coefMat,
		  const LocPointMapped& lpm ) {
  	  EXCEPTION("CoefFunctionSurfVWP:GetTensor not implemented");
  }

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<Double>& coefVec,
                 const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<Complex>& coefVec,
                 const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(Double& coefScalar,
                 const LocPointMapped& lpm ) {
  	  EXCEPTION("CoefFunctionSurfVWP:GetScalar not implemented");
  }

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(Complex& coefScalar,
		  const LocPointMapped& lpm ) {
  	  EXCEPTION("CoefFunctionSurfVWP:GetScalar not implemented");
  };

//  virtual void SetNeighborRegionId(RegionIdType id);

private:

  //! coef-function as defined in PDE
  std::map<SolutionType, shared_ptr<CoefFunctionMulti> > matCoef_;
};

} // end of namespace
#endif
