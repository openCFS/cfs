// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PLAINMATERIALHANDLER
#define FILE_PLAINMATERIALHANDLER

/**************************************************************************/
/* File:   PlainMaterialHandler.hh                                        */
/* Author: Fred Hofer                                                     */
/* Date:   14. Dez. 99                                                    */
/*                                                                        */
/* Reads Material Data from a Materialfile.                               */
/**************************************************************************/


#include "MaterialHandler.hh"

namespace CoupledField
{

  //! Class for reading material data file
  /*! 
    Includes all functions for reading material data out of file
  */


  class PlainMaterialHandler : public MaterialHandler {
  
  public:

    //! Constructor
    PlainMaterialHandler( const std::string & fileName );

    //! Loads the specified material

    //! This method loads the given material from the material file and
    //! assigns it the given pointer.
    //! \param material Empty pointer to material object
    //! \param matName Name of the material to be read
    //! \param matClass Materialclass the material belongs to
    BaseMaterial *  LoadMaterial( const std::string matName, 
                                  const MaterialClass matClass );
    
  private:
    /// pointer to a vector which holds all material data
    //    Vector<MaterialData>* mat2material;

    /// scaling with 1e-10 stiffness and 1e10 electrical part of material data
    int scaleMatDat;
        
    /// read piezolectric material data
    void ReadPiezo(std::ifstream & fin, BaseMaterial * material, bool & matC);

    /// read mechanic material data
    void ReadMechanic(std::ifstream & fin, BaseMaterial * material, 
		      bool & matC);

    /// read fluid material data
    void ReadAcoustic(std::ifstream & fin, BaseMaterial* material);

    /// read thermic material data -> e.g. for heat conduction simulation
    void ReadThermic(std::ifstream & fin, BaseMaterial * material);

    /// read thermic material data -> e.g. for heat conduction simulation
    void ReadFlow(std::ifstream & fin, BaseMaterial * material);

    /// read linear magnetic material data (could be also done by ReadMagNonLin -> is only here for compatibility)
    void ReadMagnetic(std::ifstream & fin, BaseMaterial * material);

    /// read linear electrostatic material data
    void ReadElectrostatic(std::ifstream & fin, BaseMaterial * material,  bool & matC);

    /// read next line of materialfile
    void ReadLine(std::ifstream & fin, char* buffer);

    //! Find next datarecord with materialname matName, material type matType in the 
    //!current line (if found) is written into buffer
    //! \return true, if material was found
    bool FindMat( std::ifstream &fin, const char* matName, char* buffer,
                  const char* matType );


    
  };
} // end namespace CoupledField
 
#endif
