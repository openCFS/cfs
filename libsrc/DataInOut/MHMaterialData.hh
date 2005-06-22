#ifndef MHMATERIAL_DATA
#define MHMATERIAL_DATA


#include "General/environment.hh"
#include "Utils/vector.hh"
#include "Matrix/matrix.hh"
#include "DataInOut/MHMaterialDataFile.hh"

namespace CoupledField
{

  //! Class for Material Data
  /*! 
    Interface class for getting material data
  */

  class MHMaterialData
  {
    // private:

  public:

    MHMaterialData(MaterialData * ptMaterial);

    
    
    //MHMaterialData(const MaterialData& mat);
    
    ~MHMaterialData();
    
    void updateMaterialData(Vector<Complex> & parameter, MaterialData * ptMaterial);
    // Domain * ptDomain;


    void calcParameterCurveAtElement(Vector<Complex> & parameter, 
                                     Matrix<Double> & parameterCoeff_, UInt element,UInt  N, 
                                     Integer delta, UInt pMax);

    //  private: 
    
    Vector<Complex> parameter_;
    Vector<Double> parameterR_;
    Vector<Double> parameterC_;
    Matrix<Double> parameterCoeff_;
    MaterialData * ptMaterial_;
    
    MHMaterialDataFile * ptMHMaterialDataFile;


  }; // end class MHMaterialData
}
 
#endif
