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
    XMLMaterialHandler( );

    //! Destructor
    ~XMLMaterialHandler();
    
    //! Load from a given file
    void LoadFromFile( const std::string& fileName );
    
    //! Load from a given string
    void LoadFromString( const std::string& str );
    
    //! Loads the specified material

    //! This method loads the given material from the material file and
    //! assigns it the given pointer.
    //! \param material Empty pointer to material object
    //! \param matName Name of the material to be read
    //! \param matClass Materialclass the material belongs to
    BaseMaterial *  LoadMaterial( const std::string &matName,
                                  MaterialClass matClass );
    
  private:

    //! Path to the material XML schema relative to share/xml
    static const std::string schemaFile_;

    //! XML Schema URL
    static const std::string schemaUrl_;

    //! The param node contains the material xml file content
    PtrParamNode rootNode_;

    //! =====================================================================
    //! METHODS FOR READING MATERIALS 
    //! =====================================================================

    //! Reads piezo material.
    //! \param material Material object to be filled with data
    void ReadPiezo(BaseMaterial *material, PtrParamNode pn);

    //! Reads mechanical material.
    //! \param material Material object to be filled with data
    void ReadMechanic(BaseMaterial *material, PtrParamNode pn);

    //! Reads smooth material.
    //! \param material Material object to be filled with data
    void ReadSmooth(BaseMaterial *material, PtrParamNode pn);

    //! Reads acoustic material.
    //! \param material Material object to be filled with data
    void ReadAcoustic(BaseMaterial *material, PtrParamNode pn);

    //! Reads electric material.
    //! \param material Material object to be filled with data
    void ReadElectrostatic(BaseMaterial *material, PtrParamNode pn); 

    //! Reads electro-quasistatic material.
    //! \param material Material object to be filled with data
    void ReadElecQuasistatic(BaseMaterial *material, PtrParamNode pn);

    //! Reads magnetic material.
    //! \param material Material object to be filled with data
    void ReadMagnetic(BaseMaterial *material, PtrParamNode pn);

    //! Reads thermic material.
    //! \param material Material object to be filled with data
    void ReadThermic(BaseMaterial *material, PtrParamNode pn);

    //! Reads flow material.
    //! \param material Material object to be filled with data
    void ReadFlow(BaseMaterial *material, PtrParamNode pn);

    //! Reads material for TEST PDE.
    //! \param material Material object to be filled with data
    void ReadTest(BaseMaterial *material, PtrParamNode pn);

    //! Reads thermoelastic material.
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadThermoelastic(BaseMaterial *material, PtrParamNode pn);

    //! Reads pyroelectric material.
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadPyroelectric(BaseMaterial *material, PtrParamNode pn);
    
    //! Reads magnetostrictive material
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadMagStrict(BaseMaterial *material, PtrParamNode pn);

    //! Reads electric conduction material
    //! \param material Material object to be filled with data
    //! \param matName Name of the material to be read
    void ReadElectricConduction(BaseMaterial *material, PtrParamNode pn);

    //! =====================================================================
    //! HELPER METHODS
    //! =====================================================================

    //! Read a scalar value
    bool ReadScalar( PtrParamNode ptrNode, std::string& val, 
                     std::string str1, std::string str2 );

    //! Read a scalar value and pass it to material class
    PtrCoefFct ReadScalar( PtrParamNode ptrNode, std::string str,
                           Global::ComplexPart type);

    //! Reads a linear scalar value
    bool ReadScalarLin( PtrParamNode ptrNode, std::string& val, 
                        std::string str1, std::string str2 );

    //! Read a scalar value and pass it to material class
    PtrCoefFct ReadScalarLin( PtrParamNode ptrNode, std::string str,
                              Global::ComplexPart type);

    //! Read tensor value and pass it to material class
    PtrCoefFct ReadTensor( PtrParamNode ptrNode, Global::ComplexPart type );
    
    //! Read square 3x3 tensor, allowing different definitions
    
    //! This method reads the definition of a 3x3 tensor, which can be defined
    //! either by an isotropic value (diaganal value), transversal-isotropic,
    //! orthotropic or as a full tensor representation.
    //! \param ptrNode    ParamNode above <isotropic>, <orthotropic> etc.
    //! \param mat        Pointer to material objects
    //! \param isoProp    Enum for the isotropic paramater
    //! \param orthoProp  3 Enums for the 3 distinct orthotropic scalar values
    //! \param tensorProp Enum for the final tensor property to be set
    //! \param part       Complex part (real, imag, complex)
    void ReadSquare3x3Tensor( PtrParamNode ptrNode, BaseMaterial *mat,
                              MaterialType isoPProp,
                              MaterialType* orthoProp,
                              MaterialType tensorProp,
                              Global::ComplexPart part = Global::COMPLEX);
    
    //! Read a mechanical stiffness tensor
    BaseMaterial::SymmetryType ReadStiffnessTensor(PtrParamNode ptrNode,
                                                   Global::ComplexPart type,
                                                   BaseMaterial::CoefMap &coefMap);

    //! Read a smooth stiffness tensor
    BaseMaterial::SymmetryType ReadStiffnessTensorSmooth(PtrParamNode ptrNode,
                                                          Global::ComplexPart type,
                                                          BaseMaterial::CoefMap &coefMap);

    //! Reads the bh values out of the material file
    void ReadXYValues( PtrParamNode paramNode,
                       const std::string &xName,
                       const std::string &yName,
                       Vector<Double>& xValues,
                       Vector<Double>& yValues );

    //! Read Rayleigh damping parameters and safe info into coefFunctions
    void ReadRayleighDamping(PtrParamNode paramNode, BaseMaterial *material);
    //! Read lossTanDelta damping and safe info into coeffunctions
    void ReadLossTanDeltaDamping(PtrParamNode paramNode, BaseMaterial *material);

    //! Read nonlinearity descriptor
    BaseMaterial::MatDescriptorNl ReadNonlinDescriptor(PtrParamNode paramNode,
                                                       BaseMaterial *material);

    //! Read hysteresis; central function to avoid endless duplicates
    void ReadHysteresis(BaseMaterial *material, PtrParamNode hystNode, bool isMagnetic);
    void ReadHystOperator(BaseMaterial *material, PtrParamNode operatorNode, bool setStrains, bool isMagnetic);
  };
  
}

#endif
