#include "singleDriver.hh"


namespace CoupledField 
{
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

  void calcMultiHarmMassMatrix(Double & omega);

private:

  Double  startFreq_;
  Double  stopFreq_;
  Integer numFreq_;
  Integer saveType_;
  Integer nrMultHarms_;


  }; // end of class piezoParamIdent
}



