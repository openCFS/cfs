//#ifndef FILE_
//#define FILE_PIEZO_PARAM_IDENT

#include "Driver/singleDriver.hh"
#include <PDE/SinglePDE.hh>
#include "singleDriver.hh"

namespace CoupledField 
{

  // forward class declarations
  class MaterialData;

  //! multiharmonic analysis
  class MultiHarmonicDriver : public SingleDriver
  {

  public:

    MultiHarmonicDriver(Domain * adomain,
                        UInt stepOffset = 0,
                        Double timeOffset = 0.0,
                        std::string driverTag = "anyTag",
                        Boolean isPartOfSequence = FALSE);

    //! Destructor
    virtual ~MultiHarmonicDriver();
  
    std::ifstream * allMeasuredData;
    std::ofstream * impedCurve;
    std::ofstream * piezoLog;
    std::ofstream * parLog;

    StdVector<RegionIdType> MHsubdoms_;  //!< subdomain-levels belongig to PDE
    StdVector<RegionIdType> MHsurfdoms_; //!< surface-domain-levels belongig to PDE

  
    virtual void SolveProblem();
    virtual void MHAssembleMatrices();

  protected:
    //! \param parameter - new set of piezoelectric material parameters
    void updateMaterialData(Vector<Double> & parameter, MaterialData * ptMaterial);
    // Domain * ptDomain;


    void calcParameterCurveAtElement(Vector<Complex> & parameter, 
                                     Matrix<Double> & parameterCoeff_, UInt element,Integer  N, 
                                     Integer delta, UInt pMax);

    SinglePDE * ptMyPDE_;

    void createMHMassMatrix(UInt N);

    // pointers to classes involved
    BaseSystem * ptAlgsys_;
    Assemble * ptAssemble_;

    Integer nrMultHarms_;

    UInt getNrMultHarms(){
      return nrMultHarms_;
    }
    

  private:

    Double  startFreq_;
    Double  stopFreq_;
    UInt numFreq_;
    UInt saveType_;

    Matrix<Double> parameterCoeff_;
    Vector<Complex> parameter_;




  }; // end of class piezoParamIdent
}



