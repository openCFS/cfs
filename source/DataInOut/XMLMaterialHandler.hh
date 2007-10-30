// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_XMLMATERIALHANDLER
#define FILE_XMLMATERIALHANDLER

#include "MaterialHandler.hh"

namespace CoupledField {

  //! Class for reading materials fron xml-file
  class XMLMaterialHandler : public MaterialHandler {
  
  public:

    //! Constructor
    XMLMaterialHandler( const std::string & fileName );

    //! Destructor
    ~XMLMaterialHandler();

    //! Loads the specified material

    //! This method loads the given material from the material file and
    //! assigns it the given pointer.
    //! \param material Empty pointer to material object
    //! \param matName Name of the material to be read
    //! \param matClass Materialclass the material belongs to
    BaseMaterial *  LoadMaterial( const std::string matName, 
                                  const MaterialClass matClass );
    
  private:

    //! The param node contains the material xml file content
    ParamNode* parser_;

    //! Reads piezo material.
    //! \param material Material object to be filled with data
    void ReadPiezo(BaseMaterial *material, ParamNode* pn);

    //! Reads mechanical material.
    //! \param material Material object to be filled with data
    void ReadMechanic(BaseMaterial *material, ParamNode* pn);

    //! Reads acoustic material.
    //! \param material Material object to be filled with data
    void ReadAcoustic(BaseMaterial *material, ParamNode* pn);

    //! Reads electric material.
    //! \param material Material object to be filled with data
    void ReadElectrostatic(BaseMaterial *material, ParamNode* pn); 

    //! Reads magnetic material.
    //! \param material Material object to be filled with data
    void ReadMagnetic(BaseMaterial *material, ParamNode* pn);

    //! Reads thermic material.
    //! \param material Material object to be filled with data
    void ReadThermic(BaseMaterial *material, ParamNode* pn);

    //! Reads flow material.
    //! \param material Material object to be filled with data
    void ReadFlow(BaseMaterial *material, ParamNode* pn);
    //! Reads thermoelastic material.
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadThermoelastic(BaseMaterial *material, ParamNode* pn);


    //! Reads pyroelectric material.
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadPyroelectric(BaseMaterial *material, ParamNode* pn);

  };
  
}

#endif
