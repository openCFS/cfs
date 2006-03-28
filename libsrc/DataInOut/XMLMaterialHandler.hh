#ifndef FILE_XMLMATERIALHANDLER
#define FILE_XMLMATERIALHANDLER

#include "MaterialHandler.hh"

namespace CoupledField {

  // forward class declaration
  class XMLParamHandler;

   //! Class for reading materials fron xml-file
  class XMLMaterialHandler : public MaterialHandler {
  
  public:

    //! Constructor
    XMLMaterialHandler( const std::string & fileName );

    //! Destructor
    ~XMLMaterialHandler();

    //! Get specific material object

    //! Loads the specified material into the given material object.
    //! \param material Material object to be fille with data
    //! \param matName Name of the material to be read
    //! \param 
    void GetMaterial( BaseMaterial * material, 
                      const std::string matName, 
                      const MaterialClass matClass );

  private:

    //! Pointer to xml-parser object
    XMLParamHandler * parser_;

    //! Reads piezo material.
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadPiezo(BaseMaterial *material,
                   const std::string matName);

    //! Reads mechanical material.
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadMechanic(BaseMaterial *material,
                      const std::string matName);

    //! Reads acoustic material.
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadAcoustic(BaseMaterial *material,
                           const std::string matName);

    //! Reads electric material.
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadElectrostatic(BaseMaterial *material,
                           const std::string matName);

    //! Reads magnetic material.
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadMagnetic(BaseMaterial *material,
                      const std::string matName);

    //! Reads thermic material.
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadThermic(BaseMaterial *material,
                      const std::string matName);

    //! Reads flow material.
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadFlow(BaseMaterial *material,
                      const std::string matName);

  };
  
}

#endif
