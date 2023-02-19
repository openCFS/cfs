#ifndef FILE_COEFFUNCTION_SURF_HH
#define FILE_COEFFUNCTION_SURF_HH


#include "CoefFunction.hh"
#include "CoefFunctionMulti.hh"
#include "FeBasis/FeFunctions.hh"
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

  //! Add a surface group
  virtual void AddEntity(shared_ptr<EntityList> entity);
  
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

  //! List of surface groups
  std::map< std::string, shared_ptr<EntityList> > entities_;

  //! Flag, if normal mapping should be performed
  bool mapNormal_;

  //! Number of DoFs (for vector-valued results)
  UInt numDofs_;
  
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


//! This class represents coefficient functions, which are defined just on a
//! surface and computes the force defined by virtual work principle
//! It's derived from CoefFunctionSurf
template<class FE>
class CoefFunctionSurfVWPnew : public CoefFunctionSurf {
  public:

    //! Constructor

    //! \param matCoef    Material parameter
    //! \param vacuumCoef Material parameter of vacuum
    //! \param surfInfo   Result info object for surface result

    CoefFunctionSurfVWPnew(PtrCoefFct matCoef, 
                           shared_ptr<ResultInfo> surfInfo, Grid* ptGrid,
                           shared_ptr<BaseFeFunction> feFnc);

    //! Destructor
    virtual ~CoefFunctionSurfVWPnew();

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

    //! Returns the total force summing up all nodal forces over a group
    void GetTotalForce(const std::string & entityName,
                       Vector<Double> & totalForce);

  private:

    //! Calculate element force
    //! \param force          (output) Array containing nodal forces
    //!                                (dim x nodes) of each element
    //! \param ptElement      (input)  Pointer to element
    //! \param dim            (input)  number of dofs = dim
    //! \param isBoundaryNode (input)  contains 1, if corresponding node is a
    //!                                boundary node, otherwise 0
    void CalcElemForce(Matrix<Double>& force, const Elem * ptElement,
                       const std::vector<bool> & isBoundaryNode);

    //! Update the cache with nodal forces
    void UpdateCache();

    //! Calculates the expression \f[ \frac{\delta \vert J \vert}{\delta r} /f]
    //! \param J    (input) Jacobian matrix
    //! \param J_dr (input) derivative of Jacobian matrix in r-direction
    //! \param lpm  LocPointMapped at which the expression shall be calculated
    static Double CalcDetJDr(const Matrix<Double> &J,
                             const Matrix<Double> &dJ_dr,
                             const LocPointMapped &lpm);

    //! pointer to the grid
    Grid* ptGrid_;

    //! CoefFunction of material parameter
    PtrCoefFct matCoef_;

    //!FeFunction
    shared_ptr<BaseFeFunction> FeFunction_;

    //! Cache for nodal forces (maps node number to force vector)
    boost::unordered_map< UInt, Vector<Double> > nodalForces_;

    //! Cache for total force per group
    std::map< std::string, Vector<Double> > totalForces_;

    //! Step number of cached data
    UInt cacheStep_;
};

} // end of namespace
#endif
