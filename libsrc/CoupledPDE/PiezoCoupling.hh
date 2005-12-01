#ifndef FILE_PIEZOCOUPLING_HH
#define FILE_PIEZOCOUPLING_HH

#include "BasePairCoupling.hh"
#include "Utils/elemstoresol.hh"
#include "Utils/nodestoresol.hh"
#include "PDE/SinglePDE.hh"



namespace CoupledField
{

  // Forward declarations
  class BaseForm;
  class MaterialData;
  class SinglePDE;

  //! Implements the definition of pairwise piezo-coupling
  class PiezoCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to first coupling PDE
    //! \param pde2 pointer to second coupling PDE
    PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2 );

    //! Destructor
    virtual ~PiezoCoupling();

    //! Trigger calculation of postprocessing results
    void PostProcess();

    void calcMaterialMatrices(Matrix<Double> &sMat, 
                              Matrix<Double>&cMat,
                              Matrix<Double> &pMat,
                              Matrix<Double> *matDat);

    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    void WriteResultsInFile(const UInt kstep,
                             const Double asteptime,
                             UInt stepOffset,
                             Double timeOffset);    
    
  protected:

    //! Obtain information on desired output quantities from parameter file
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    void ReadStoreResults();
    
    //! Definition of the (bi)linear forms
    void DefineIntegrators();
    
    //! Get correct stiffness integrator
    BaseForm * GetStiffIntegrator( MaterialData * actSDMat,
                                   Boolean reducedInt = FALSE ,
                                   Boolean isdamping = FALSE );

    // Data section
    Boolean hasOutput_;

  private:
    
    AnalysisType analysistype_;

    //! Postprocession section

    //! Contains the subdomains, on which the electric field is computed
    StdVector<RegionIdType> calcDfield_;

    //! conatins electric flux density
    BaseElemStoreSol *  Dfield_;  

    //! contains the subdomains, on which the stress is computed
    StdVector<RegionIdType> calcStress_;  

    //! conatins mechanical stresses
    BaseElemStoreSol * stress_;

    //! contains electric charges
    BaseElemStoreSol * charges_;

    //! contains neighbouring regions for charge calculation
    StdVector<RegionIdType> chargeNeighborRegion_;
    StdVector<RegionIdType> calcCharge_;

    //! contains mechanic velocity
    NodeStoreSol<Double> solDeriv1_;
    
    //! contains mechanic acceleration
    NodeStoreSol<Double> solDeriv2_;

    //! computes stresses, i.e. \sigma = cBu + e \grad \phi
    template <class TYPE>
    void CalcStress();

    //! calculate Charges
    template <class TYPE>
    void CalcCharges();

    //!! calculate electric field
    //    void CalcDfield();

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class PiezoCoupling
  //! 
  //! \purpose 
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace

#endif
