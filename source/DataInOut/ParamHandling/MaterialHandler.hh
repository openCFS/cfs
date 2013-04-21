// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MATERIALHANDLER_HH
#define FILE_MATERIALHANDLER_HH


#include "General/Environment.hh"
#include "Materials/BaseMaterial.hh"
#include "Domain/Domain.hh"
#define bufLength 200

namespace CoupledField
{

  //! Base class for reading material information
  class MaterialHandler {
 
  public:

    //! Constructor
    MaterialHandler( const std::string & fileName ) {
      fileName_ = fileName;
      domain_ = NULL;
      mp_ = NULL;
    }
    
    //! Destructor
    virtual ~MaterialHandler() {};
    
    //! Set domain
    void SetDomain(Domain* domain ) {
      this->domain_ = domain;
      this->mp_ = domain_->GetMathParser();
    }
    
    //! Retrieve name of material file
    std::string GetFileName() const {
      return fileName_;
    }
    
    //! Loads the specified material
    
    //! This method loads the given material from the material file and
    //! assigns it the given pointer.
    //! \param material Empty pointer to material object
    //! \param matName Name of the material to be read
    //! \param matClass Materialclass the material belongs to
    virtual BaseMaterial *  LoadMaterial( const std::string matName, 
                                          const MaterialClass matClass ) = 0;

  protected:
    
    //! Private standard constructor
    MaterialHandler() {};

    //! Filename of material file
    std::string fileName_;
    
    //! Pointer to domain
    Domain* domain_;
    
    //! Pointer to MathParser
    MathParser* mp_;

  };



  
#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MaterialHandler
  //! 
  //! \purpose 
  //! This class defines an abstract interface for reading material information
  //! from external sources into own material representation within CFS
  //! 
  //! \collab 
  //! Uses the interface defined by BaseMaterial
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif


} // end namespace CoupledField
#endif
