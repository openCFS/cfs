#ifndef FILE_LOADMATERIALDATA
#define FILE_LOADMATERIALDATA

/**************************************************************************/
/* File:   LoadMaterialData.hh                                            */
/* Author: Fred Hofer                                                     */
/* Date:   14. Dez. 99                                                    */
/*                                                                        */
/* Reads Material Data from a Materialfile.                               */
/**************************************************************************/


#include <Utils/vector.hh>
#include <General/environment.hh>

#include "MaterialData.hh"


#define bufLength 200

namespace CoupledField
{

  //! Class for reading material data file
  /*! 
    Includes all functions for reading material data out of file
  */


  class LoadMaterialData
  {
  public:

    /// constructor: loads all materials defined in the parserfile from the materialfile "filename"
    //    LoadMaterialData(const std::string fileName);
    LoadMaterialData(const char * fileName){};

    /// private constructor, loads all materials defined in the parserfile from the materialfile "filename"
    virtual void GetMaterial( MaterialData& material, const std::string matName, const std::string matType)=0;
    
    /// rotate into decired euler angles
    //void LoadMaterialData :: EulerAnglesRotate(MaterialData * material, const Vector<Double> & eulerAngles)

    // read nonlinear magnetic record
    // void ReadMagNonLin(std::ifstream & fin, MaterialData * material);

  };
} // end namespace CoupledField
 
#endif
