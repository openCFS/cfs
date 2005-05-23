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
                        Integer stepOffset = 0,
                        Double timeOffset = 0.0,
                        std::string driverTag = "anyTag",
                        Boolean isPartOfSequence = FALSE);

    //! Destructor
    virtual ~MultiHarmonicDriver();
  
    std::ifstream * allMeasuredData;
    std::ofstream * impedCurve;
    std::ofstream * piezoLog;
    std::ofstream * parLog;
  
    virtual void SolveProblem();

  protected:
    //! \param parameter - new set of piezoelectric material parameters
    void updateMaterialData(Vector<Double> & parameter, MaterialData * ptMaterial);
    // Domain * ptDomain;

    SinglePDE * ptMyPDE_;

    void createMHMassMatrix(Integer N);

    // pointers to classes involved
    //  SinglePDE * ptMyPDE_;
    BaseSystem * ptAlgsys_;
    Assemble * ptAssemble_;


  private:

    Double  startFreq_;
    Double  stopFreq_;
    Integer numFreq_;
    Integer saveType_;
    Integer nrMultHarms_;




  }; // end of class piezoParamIdent
}



