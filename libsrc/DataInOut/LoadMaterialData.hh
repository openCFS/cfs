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
  private:
    /// pointer to a vector which holds all material data
    //    Vector<MaterialData>* mat2material;

    /// name of materialfile
    const char * filename;    

    /// Euler Angles for rotation of anisotropic matrials
    //    Vector< Vector<Double>* > eulerAngles;


    /// scaling with 1e-10 stiffness and 1e10 electrical part of material data
    int scaleMatDat;
	
    // read nonlinear magnetic record
    // void ReadMagNonLin(std::ifstream & fin, MaterialData * material);

    /// read piezolectric material data
    void ReadPiezo(std::ifstream & fin, MaterialData * material);

    /// read fluid material data
    void ReadFluid(std::ifstream & fin, MaterialData * material);

    /// read linear magnetic material data (could be also done by ReadMagNonLin -> is only here for compatibility)
    void ReadMagnetic(std::ifstream & fin, MaterialData * material);

    /// read next line of materialfile
    void ReadLine(std::ifstream & fin, char* buffer);

    /// find next datarecord with materialname matName, the current line (if found) is written into buffer
    //    void FindMat(std::ifstream & fin, const char* matName, char* buffer);
    void FindMat(std::ifstream & fin, const char* matName, char* buffer);

    /// rotate into decired euler angles
    void EulerAnglesRotate(MaterialData * material, const Vector<Double> & eulerAngles);

  public:

    /// constructor: loads all materials defined in the parserfile from the materialfile "filename"
    //    LoadMaterialData(const std::string fileName);
    LoadMaterialData(const char * fileName);


    /// private constructor, loads all materials defined in the parserfile from the materialfile "filename"
    void GetMaterial( MaterialData& material, const std::string matName, const std::string matType);
    
  };
} // end namespace CoupledField
 
#endif
