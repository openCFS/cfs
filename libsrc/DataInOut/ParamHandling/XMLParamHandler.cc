#ifdef USE_XERCES

// includes for Xerces
#include "xercesc/util/PlatformUtils.hpp"
#include "xercesc/dom/DOM.hpp"
#include "xercesc/dom/DOMImplementation.hpp"
#include "xercesc/dom/DOMImplementationLS.hpp"
#include "xercesc/dom/DOMWriter.hpp"
#include "xercesc/framework/StdOutFormatTarget.hpp"
#include "xercesc/framework/LocalFileFormatTarget.hpp"
#include "xercesc/parsers/XercesDOMParser.hpp"
#include "xercesc/util/XMLUni.hpp"

// we want to use the Xerces-C++ namespace
using namespace xercesc;

// we want to use new with nothrow
#include <new>

// we need some character conversion routines
#include <string>

// we want to use classical C IO-routines
#include <stdio.h>
#include <fstream>

// includes from CFS++ itself
#include "General/environment.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/XMLParamHandler.hh"
#include "DataInOut/ParamHandling/XMLParserErrorHandler.hh"
#include "Utils/tools.hh"


// uncomment for debugging (shows details of search progress on screen)
// this does not use (*debug), so it works also with a non-debug version
// #define DEBUG_XMLPARAMHANDLER


namespace CoupledField {


  // *************************************************************************
  //   Constructors and Destructors
  // *************************************************************************


  // =========================================
  //   Constructor for XMLParamHandler class
  // =========================================
  XMLParamHandler::XMLParamHandler( const char *fname ): parser_(NULL) {

    ENTER_FCN( "XMLParamHandler::XMLParamHandler" );

    // Initialise the XML4C2 system
    try {
      XMLPlatformUtils::Initialize();
    }
    catch( const XMLException &event ) {
      (*error) << "The following error occured during initialisation of "
               << "xerces-c!\n"
               << event.getMessage();
      Error( __FILE__, __LINE__ );
    }


#ifdef DEBUG
    std::cout << "\n\n XML parsers uses the Schema: http://www.cfs++.org\n "
              << commandLine->GetSchemaPath() << "/CFS.xsd" << std::endl;
#endif


    // Generate progress report (tell user what we are doing)
    std::stringstream msg;
    msg << "Reading parameters from file '" << fname << "'";
    Info->StartProgress( msg.str() );

    // Check if file exists
    std::ifstream inFile;
    inFile.open( fname );
    if ( inFile == NULL ) {
      (*error) << "The file '" << fname << "' could not be opened!";
      Error( __FILE__, __LINE__ );
    }
    inFile.close();
 
    // Generate parser and parse XML parameter file
    rootElem_ = ParseFile( &parser_, fname );

    // Generate parser and parse XML defaults file
    cfsDefaults_ = commandLine->GetSchemaPath();
    cfsDefaults_ += "/Defaults/CFS++Defaults.xml";
    rootElemDefaults_ = ParseFile( &parserDefaults_, cfsDefaults_.c_str() );

    // Toggle verbosity
#ifdef DEBUG_XMLPARAMHANDLER
    beVerbose_ = true;
#else
    beVerbose_ = false;
#endif

    // That's it
    Info->FinishProgress();

  }


  // ========================================
  //   Destructor for XMLParamHandler class
  // ========================================
  XMLParamHandler::~XMLParamHandler() {

    ENTER_FCN( "XMLParamHandler::~XMLParamHandler" );

#ifdef DEBUG_XMLPARAMHANDLER
    PrintTree( parser_->getDocument() );
    PrintTree( parserDefaults_->getDocument() );
#endif

    // Delete parser for parameter file
    if ( parser_ != NULL ) {
      delete parser_->getErrorHandler();
      delete parser_;
    }

    // Delete parser for defaults' file
    if ( parserDefaults_ != NULL ) {
      delete parserDefaults_->getErrorHandler();
      delete parserDefaults_;
    }

    // Shutdown platform dependend utilities
    XMLPlatformUtils::Terminate();

  }


  // =================================================
  //   Default constructor for XMLParamHandler class
  // =================================================
  XMLParamHandler::XMLParamHandler() {}


  // *************************************************************************
  //   Public Methods: Simple Query Functions
  // *************************************************************************


  // ================================================================
  //   Return as string the value of a parameter matching a keyword
  // ================================================================
  void XMLParamHandler::Get( const std::string key, std::string &value,
                             const std::string section,
                             const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // Generate vectors of keywords and side-constraints
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    GenerateSearchParams( key, section, subsection, keyVec, attrVec, valVec );

    // Use constrained get to find the value
    Get( keyVec, attrVec, valVec, value );

    // Tell what we found
    if ( beVerbose_ == true ) {
      std::string msg = "Get: Value for parameter '" + key;
      msg += "' is '" + value + "'";
      Info->Warning( msg );
    }

  }


  // =================================================================
  //   Return as integer the value of a parameter matching a keyword
  // =================================================================
  void XMLParamHandler::Get( const std::string key, int &value,
                             const std::string section,
                             const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // Find element/value matching keyword as string
    std::string match;
    Get( key, match, section, subsection );

    // Error detection and default handling occured in the above Get,
    // so we only need to convert the detected value
    value = atoi( match.c_str() );

    // Tell what we found
    if ( beVerbose_ == true ) {
      std::string msg = "Get: Value for parameter '" + key;
      msg += "' is '" + match + "'";
      Info->Warning( msg );
    }

  }


  // =================================================================
  //   Return as integer the value of a parameter matching a keyword
  // =================================================================
  void XMLParamHandler::Get( const std::string key, UInt &value,
                             const std::string section,
                             const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // Find element/value matching keyword as string
    std::string match;
    Get( key, match, section, subsection );

    // Error detection and default handling occured in the above Get,
    // so we only need to convert the detected value
    value = String2UInt( match.c_str() );

    // Tell what we found
    if ( beVerbose_ == true ) {
      std::string msg = "Get: Value for parameter '" + key;
      msg += "' is '" + match + "'";
      Info->Warning( msg );
    }

  }


  // ================================================================
  //   Return as double the value of a parameter matching a keyword
  // ================================================================
  void XMLParamHandler::Get( const std::string key, double &value,
                             const std::string section,
                             const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // Find element/value matching keyword as string
    std::string match;
    Get( key, match, section, subsection );

    // Error detection and default handling occured in the above Get,
    // so we only need to convert the detected value
    value = atof( match.c_str() );

    // Tell what we found
    if ( beVerbose_ == true ) {
      std::string msg = "Get: Value for parameter '" + key;
      msg += "' is '" + match + "'";
      Info->Warning( msg );
    }

  }


  // ===================================================
  //   Return list of string values matching a keyword
  // ===================================================
  void XMLParamHandler::GetList( const std::string key,
                                 StdVector<std::string> &list,
                                 const std::string section,
                                 const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::GetList" );

    // Clear input vector
    ClearVector( list );

    // Generate vectors of keywords and side-constraints
    StdVector<std::string> keys;
    StdVector<std::string> attribs;
    StdVector<std::string> values;
    GenerateSearchParams( key, section, subsection, keys, attribs, values );

    // Now we just have to pass everything through to the
    // central auxilliary search method
    FindAllMatches( keys, attribs, values, list, rootElem_ );

  }


  // ===================================================
  //   Return list of Double values matching a keyword
  // ===================================================
  void XMLParamHandler::GetList( const std::string key,
                                 StdVector<Double> &list,
                                 const std::string section,
                                 const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::GetList" );

    // Clear input vector
    ClearVector( list );


    // First determine all matches as strings
    StdVector<std::string> matches;
    GetList( key, matches, section, subsection );

    // Now perform conversion
    for ( unsigned int i = 0; i < matches.GetSize(); i++ ) {
      list.Push_back( String2Double( matches[i] ) );
    }
  }


  // ====================================================
  //   Return list of Integer values matching a keyword
  // ====================================================
  void XMLParamHandler::GetList( const std::string key,
                                 StdVector<Integer> &list,
                                 const std::string section,
                                 const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::GetList" );

    // Clear input vector
    ClearVector( list );


    // First determine all matches as strings
    StdVector<std::string> matches;
    GetList( key, matches, section, subsection );

    // Now perform conversion
    for ( unsigned int i = 0; i < matches.GetSize(); i++ ) {
      list.Push_back( String2Integer( matches[i] ) );
    }
  }


  // =================================================
  //   Return list of UInt values matching a keyword
  // =================================================
  void XMLParamHandler::GetList( const std::string key,
                                 StdVector<UInt> &list,
                                 const std::string section,
                                 const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::GetList" );

    // Clear input vector
    ClearVector( list );


    // First determine all matches as strings
    StdVector<std::string> matches;
    GetList( key, matches, section, subsection );

    // Now perform conversion
    for ( unsigned int i = 0; i < matches.GetSize(); i++ ) {
      list.Push_back( String2UInt( matches[i] ) );
    }
  }


  // ================================================================
  //   Return as string the value of a parameter matching a keyword
  // ================================================================
  void XMLParamHandler::Get( const StdVector<std::string> &keyVec,
                             std::string &value ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // Generate vectors for the side-constraints
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    for ( unsigned int k = 0; k < keyVec.GetSize() - 1; k++ ) {
      attrVec.Push_back( "" );
      valVec.Push_back( "" );
    }

    // Use constrained Get to find the value
    Get( keyVec, attrVec, valVec, value );
  }


  // ================================================================
  //   Return as Double the value of a parameter matching a keyword
  // ================================================================
  void XMLParamHandler::Get( const StdVector<std::string> &keyVec,
                             Double &value ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // First get value as string
    std::string match;
    Get( keyVec, match );

    // Error detection and default handling occured in the above Get,
    // so we only need to convert the detected value
    value = atof( match.c_str() );
  }


  // =================================================================
  //   Return as Integer the value of a parameter matching a keyword
  // =================================================================
  void XMLParamHandler::Get( const StdVector<std::string> &keyVec,
                             Integer &value ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // First get value as string
    std::string match;
    Get( keyVec, match );

    // Error detection and default handling occured in the above Get,
    // so we only need to convert the detected value
    value = atoi( match.c_str() );
  }


  // ==============================================================
  //   Return as UInt the value of a parameter matching a keyword
  // ==============================================================
  void XMLParamHandler::Get( const StdVector<std::string> &keyVec,
                             UInt &value ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // First get value as string
    std::string match;
    Get( keyVec, match );

    // Error detection and default handling occured in the above Get,
    // so we only need to convert the detected value
    value = String2UInt( match.c_str() );
  }


  // *************************************************************************
  //   Public Methods: Constrained Query Functions
  // *************************************************************************


  // ================================================================
  //   Return as string the value of a parameter matching a keyword
  // ================================================================
  void XMLParamHandler::Get( const StdVector<std::string> &keyVec,
                             const StdVector<std::string> &attrVec,
                             const StdVector<std::string> &valVec,
                             std::string &value ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // Find all elements/values matching keyword in tree
    StdVector<std::string> matches;
    GetList( keyVec, attrVec, valVec, matches );

    // If there was no unique match, call problem handler
    if ( matches.GetSize() > 1 ) {
      MultipleMatchHandler( keyVec, attrVec, valVec, matches.GetSize() );
    }

    // If there was no match at all, call problem handler
    else if ( matches.GetSize() == 0 ) {
      NoMatchHandler( keyVec, attrVec, valVec, value );
    }
    
    // There was a unique match, so convert detected value
    else {
      value = matches[0];
    }
  }


  // =================================================================
  //   Return as integer the value of a parameter matching a keyword
  // =================================================================
  void XMLParamHandler::Get( const StdVector<std::string> &keyVec,
                             const StdVector<std::string> &attrVec,
                             const StdVector<std::string> &valVec,
                             Integer &value ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // Find element/value matching keyword as string
    std::string match;
    Get( keyVec, attrVec, valVec, match );

    // Error detection and default handling occured in the above Get,
    // so we only need to convert the detected value
    value = atoi( match.c_str() );
  }


  // ==============================================================
  //   Return as UInt the value of a parameter matching a keyword
  // ==============================================================
  void XMLParamHandler::Get( const StdVector<std::string> &keyVec,
                             const StdVector<std::string> &attrVec,
                             const StdVector<std::string> &valVec,
                             UInt &value ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // Find element/value matching keyword as string
    std::string match;
    Get( keyVec, attrVec, valVec, match );

    // Error detection and default handling occured in the above Get,
    // so we only need to convert the detected value
    value = (UInt) atoi( match.c_str() );
  }


  // ================================================================
  //   Return as double the value of a parameter matching a keyword
  // ================================================================
  void XMLParamHandler::Get( const StdVector<std::string> &keyVec,
                             const StdVector<std::string> &attrVec,
                             const StdVector<std::string> &valVec,
                             Double &value ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // Find element/value matching keyword as string
    std::string match;
    Get( keyVec, attrVec, valVec, match );

    // Error detection and default handling occured in the above Get,
    // so we only need to convert the detected value
    value = atof( match.c_str() );
  }


  // ===================================================
  //   Return list of string values matching a keyword
  // ===================================================
  void XMLParamHandler::GetList(  const StdVector<std::string> &keyVec,
                                  const StdVector<std::string> &attrVec,
                                  const StdVector<std::string> &valVec,
                                  StdVector<std::string> &list ) {

    ENTER_FCN( "XMLParamHandler::GetList" );

    // Clear input vector
    ClearVector( list );

    // We just have to pass everything through to the
    // central auxilliary search method
    FindAllMatches( keyVec, attrVec, valVec, list, rootElem_ );

  }


  // ===================================================
  //   Return list of Double values matching a keyword
  // ===================================================
  void XMLParamHandler::GetList(  const StdVector<std::string> &keyVec,
                                  const StdVector<std::string> &attrVec,
                                  const StdVector<std::string> &valVec,
                                  StdVector<Double> &list ) {

    ENTER_FCN( "XMLParamHandler::GetList" );

    // Clear input vector
    ClearVector( list );


    // First determine all matches as strings
    StdVector<std::string> matches;
    FindAllMatches( keyVec, attrVec, valVec, matches, rootElem_ );

    // Now perform conversion
    for ( unsigned int i = 0; i < matches.GetSize(); i++ ) {
      list.Push_back( String2Double( matches[i] ) );
    }
  }


  // ====================================================
  //   Return list of Integer values matching a keyword
  // ====================================================
  void XMLParamHandler::GetList(  const StdVector<std::string> &keyVec,
                                  const StdVector<std::string> &attrVec,
                                  const StdVector<std::string> &valVec,
                                  StdVector<Integer> &list ) {

    ENTER_FCN( "XMLParamHandler::GetList" );

    // Clear input vector
    ClearVector( list );


    // First determine all matches as strings
    StdVector<std::string> matches;
    FindAllMatches( keyVec, attrVec, valVec, matches, rootElem_ );

    // Now perform conversion
    for ( unsigned int i = 0; i < matches.GetSize(); i++ ) {
      list.Push_back( String2Integer( matches[i] ) );
    }
  }


  // =============================================================
  //   Return list of unsigned Integer values matching a keyword
  // =============================================================
  void XMLParamHandler::GetList(  const StdVector<std::string> &keyVec,
                                  const StdVector<std::string> &attrVec,
                                  const StdVector<std::string> &valVec,
                                  StdVector<UInt> &list ) {

    ENTER_FCN( "XMLParamHandler::GetList" );

    // Clear input vector
    ClearVector( list );


    // First determine all matches as strings
    StdVector<std::string> matches;
    FindAllMatches( keyVec, attrVec, valVec, matches, rootElem_ );

    // Now perform conversion
    for ( unsigned int i = 0; i < matches.GetSize(); i++ ) {
      list.Push_back( String2UInt( matches[i] ) );
    }
  }


  // *************************************************************************
  //   Public Methods: Specialised Query Functions
  // *************************************************************************


  // =====================================
  //   Return a list of the defined PDEs
  // =====================================
  void XMLParamHandler::GetPDEList( StdVector<std::string> &list ) {

    ENTER_FCN( "XMLParamHandler::GetPDEList" );

    // string for assembling error messages
    std::string errmsg;

    // Check if vector is empty. If not issue a warning
    // and erase its entries, if this is desired
    ClearVector( list );

    // Find PDE section
    XMLCh *tmpString = C2X( "pdeList" );
    DOMNodeList *pdesec = rootElem_->getElementsByTagName( tmpString );
    FreeX( &tmpString );

    // Check that there is only one such section
    if ( pdesec->getLength() != 1 ) {
      (*error) << "Got " << pdesec->getLength() << " pdeList elements "
               << "in parameter file!";
      Error( __FILE__, __LINE__ );
    }

    // The names of the PDEs are the tags of the child elements of the
    // PDE_list element, so get hold of them and make sure, there is
    // at least one child/PDE
    DOMNodeList *pdelist = pdesec->item(0)->getChildNodes();
    if ( pdelist->getLength() == 0 ) {
      (*error) << "Cannot find a single PDE in parameter file!";
      Error( __FILE__, __LINE__ );
    }

    // Now get hold of tags, convert them to strings and assemble vector
    std::string pdename;
    char *auxString = NULL;
    for ( unsigned int i = 0; i < pdelist->getLength(); i++ ) {

      // Only treat element children and not comments!
      if ( pdelist->item(i)->getNodeType() == DOMNode::ELEMENT_NODE ) {
        auxString = X2C( Node2Elem( pdelist->item(i) )->getNodeName() );
        pdename.assign( auxString );
        FreeC( &auxString );
        list.Push_back( pdename );
      }
    }
  }


  // ===========================================
  //   Return a list of iterative coupled PDEs
  // ===========================================
  void XMLParamHandler::GetIterCoupledPDEList( StdVector<std::string> &list,
                                               const std::string sequenceTag){

    ENTER_FCN( "XMLParamHandler::GetIterCoupledPDEList" );

    // Check if vector is empty. If not issue a warning
    // and erase its entries, if this is desired
    ClearVector( list );

    // Get all coupling sections in the param file
    XMLCh *tmpString = C2X( "couplingList" );
    DOMNodeList *coupledSections = rootElem_->getElementsByTagName(tmpString);
    FreeX( &tmpString );

    // Pick that coupling section, which matches
    // the specfifed sequenceTag
    DOMElement *auxElem = NULL;
    DOMElement *currentCouplingSec = NULL;
    Boolean sectionFound = FALSE;

    for (unsigned int i=0; i<coupledSections->getLength(); i++) {      
      auxElem = Node2Elem( coupledSections->item(i) );
      if (AttribHasValue( auxElem, "tag", sequenceTag, false) ) {
        // Ensure that only one section matches
        if ( sectionFound == FALSE ) {
          sectionFound = TRUE;
          currentCouplingSec = auxElem;
        }
        else {
          (*error) << "Got more than one matching coupling section for "
                   << "tag '" << sequenceTag
                   << "'.\n Please correct parameter file!";
          Error( __FILE__, __LINE__ );
        }
      }
    }

    // Print error if specified coupling section was not found
    if ( sectionFound == FALSE ) {
      (*error) << "The coupling section with tag '" << sequenceTag
               << "' was not found in the parameter file!";
      Error( __FILE__, __LINE__ );
    }

    // Get the iterative coupling section in the current
    // section of couplings
    tmpString = C2X( "iterative" );
    DOMNodeList *coupledPDEsec =
      currentCouplingSec->getElementsByTagName( tmpString );
    FreeX( &tmpString );

    // If no such section exists, simply return an empty vector
    if ( coupledPDEsec->getLength() == 0 ) {
      return;
    }

    // Check that there is only one such section
    if ( coupledPDEsec->getLength() > 1 ) {
      (*error) << "Got " << coupledPDEsec->getLength()
               << " couplingList elements in parameter file!";
      Error( __FILE__, __LINE__ );
    }

    // Find iterative Coupled section
    DOMNodeList * iterCoupledPDEsec = coupledPDEsec->item(0)->getChildNodes();
    if ( iterCoupledPDEsec->getLength() == 0 ) {
      (*error) << "Cannot find an iterative coupling section in parameter "
               << "file!";
      Error( __FILE__, __LINE__ );
    }

    // iterate over all pairwise couplings
    std::string *pdename = NULL;
    for ( unsigned int i = 0; i < iterCoupledPDEsec->getLength(); i++ ) {

      // Only treat element children and not comments!
      if ( iterCoupledPDEsec->item(i)->getNodeType() ==
           DOMNode::ELEMENT_NODE ) {
    
        // The names of the PDEs are the tags of the child elements of the
        // PDE_list element, except perhaps the last one, which specifies
        // the nonlinear coupling
        DOMNodeList *iterPDElist =
          iterCoupledPDEsec->item(i)->getChildNodes();
        if ( iterPDElist->getLength() == 0 ) {
          (*error) << "Cannot find a single PDE in iterative couplingList "
                   << " of parameter file!";
          Error( __FILE__, __LINE__ );
        }

        // Now get hold of tags, convert them to strings and assemble vector
        Boolean found = FALSE;
        for ( unsigned int k = 0; k < iterPDElist->getLength(); k++ ) {
 
          // Only treat element children and not comments!
          if ( iterPDElist->item(k)->getNodeType() == DOMNode::ELEMENT_NODE ){

            pdename = X2S( Node2Elem( iterPDElist->item(k) )->getNodeName() );

            // only get elements which describe a PDE element
            if ( *pdename != "nonLinear" ) {

              // Now ensure, that each PDEname occurs only one time
              found = FALSE;
              for ( unsigned int j = 0; j < list.GetSize(); j++ ) {
                if ( *pdename == list[j] ) {
                  found = TRUE;
                  break;
                }
              }
              if ( !found ) {
                list.Push_back( *pdename );
              }
            }

            // Clean up
            FreeS( &pdename );
          }
        }
      }
    }
  }
  

  // =========================
  //   GetDirectCouplingList
  // =========================
  void
  XMLParamHandler::GetDirectCouplingList( StdVector<std::string> &list,
                                          const std::string sequenceTag ) {

    ENTER_FCN( "XMLParamHandler::GetDirectCouplingList" );

    // Check if vector is empty. If not issue a warning
    // and erase its entries, if this is desired
    ClearVector( list );

    // Get all coupling sections in the param file
    XMLCh *tmpString = C2X( "couplingList" );
    DOMNodeList *coupledSections = rootElem_->getElementsByTagName(tmpString);
    FreeX( &tmpString );

    // Pick that coupling section, which matches
    // the specfifed sequenceTag
    DOMElement *auxElem = NULL;
    DOMElement *currentCouplingSec = NULL;
    Boolean sectionFound = FALSE;
    
    for ( unsigned int i = 0; i < coupledSections->getLength(); i++) {
      auxElem = Node2Elem( coupledSections->item(i) );
      if (AttribHasValue( auxElem, "tag", sequenceTag, false) ) {
        // Ensure that only one section matches
        if ( sectionFound == FALSE ) {
          sectionFound = TRUE;
          currentCouplingSec = auxElem;
        }
        else {
          (*error) << "Got more than one matching coupling section for tag '"
                   << sequenceTag << "'.\n Please correct parameter file!";
          Error( __FILE__, __LINE__ );
        }
      }
    }
     
    // Print error if specified coupling section was not found
    if ( sectionFound == FALSE ) {
      (*error) << "The coupling section with tag '" << sequenceTag
               << "' was not found in the parameter file!";
      Error( __FILE__, __LINE__ );
    }

    // Get the iterative coupling section in the current
    // section of couplings
    tmpString = C2X( "direct" );
    DOMNodeList *coupledPDEsec =
      currentCouplingSec->getElementsByTagName( tmpString );
    FreeX( &tmpString );

    // If no direct coupled section is found, return simply an 
    // empty vector

    if ( coupledPDEsec->getLength() == 0 ) {
      return;
    }

    // Check that there is at most one such section
    if ( coupledPDEsec->getLength() > 1 ) {
      (*error) << "Got " << coupledPDEsec->getLength()
               << " 'direct' elements in parameter file!\n"
               << "One at maximum is allowed in the 'coupledList' section!";
      Error( __FILE__, __LINE__ );
    }

    // Find direct Coupled section
    DOMNodeList *directCoupledPDEsec =
      coupledPDEsec->item(0)->getChildNodes();

    if ( directCoupledPDEsec->getLength() == 0 ) {
      (*error) << "Cannot find a direct coupling section in parameter file!";
      Error( __FILE__, __LINE__ );
    }

    // iterate over all pairwise couplings
    std::string pdename;
    char *auxString = NULL;
    for ( unsigned int i = 0; i < directCoupledPDEsec->getLength(); i++ ) {

      // Only treat element children and not comments!
      if ( directCoupledPDEsec->item(i)->getNodeType() ==
           DOMNode::ELEMENT_NODE){

        // The names of the PDEs are the tags of the child elements of the
        // PDE_list element, except perhaps the last one, which specifies
        // the nonlinear coupling
        DOMNodeList *directCouplinglist = 
          directCoupledPDEsec->item(i)->getChildNodes();

        if ( directCouplinglist->getLength() == 0 ) {
          (*error) << "Cannot find a single direct coupling definition in "
                   << "the direct coupling section of the parameter file!";
          Error( __FILE__, __LINE__ );
        }
  
        // Now get hold of name, convert it and push it back in the name list
        auxString =
          X2C( Node2Elem( directCoupledPDEsec->item(i) )->getNodeName() );
        pdename.assign( auxString );
        list.Push_back( pdename );
        FreeC( &auxString );
      }
    }
  }

    
  // ======================================
  //   Return a list of the defined coils
  // ======================================
  void XMLParamHandler::GetCoilList( StdVector<std::string> &list,
                                     const std::string pde,
                                     const std::string sequenceTag ) {

    ENTER_FCN( "XMLParamHandler::GetCoilList" );

    // Check if vector is empty. If not issue a warning
    // and erase its entries, if this is desired
    ClearVector( list );

    // Assemble keywords for attribute search
    // and empty vectors for side constraints
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    keyVec.Push_back( "pdeList" );
    attrVec.Push_back( "" );
    valVec.Push_back( "" );
    if ( pde != "" ) {
      keyVec.Push_back( pde );
      attrVec.Push_back( "" );
      valVec.Push_back( "" );
    }

    keyVec.Push_back( "coils" );
    attrVec.Push_back( "" );
    valVec.Push_back( "" );

    keyVec.Push_back( "*" );
    attrVec.Push_back( "tag" );
    valVec.Push_back( sequenceTag );
    
    // Find coil names
    StdVector<DOMAttr*> *attrs =
      FindMatchingAttributes( "name", keyVec, attrVec, valVec, rootElem_ );

    char *value = NULL;
    for ( unsigned int i = 0; i < attrs->GetSize(); i++ ) {
      value = X2C( (*attrs)[i]->getValue() );
      list.Push_back( value );
      FreeC( &value );
    }

    // Cleanup
    if ( attrs != NULL ) {
      attrs->Clear();
      delete attrs;
    }
  }


  // =====================================
  //   Return the type of a certain coil
  // =====================================
  void XMLParamHandler::GetCoilType( std::string &coilType,
                                     const std::string coilName,
                                     const std::string pde ) {

    ENTER_FCN( "XMLParamHandler::GetCoilType" );

    // Be verbose if desired
    if ( beVerbose_ == true ) {
      std::cerr << "\n Trying to get description of coil named '"
                << coilName << "'\n\n";
    }

    // First generate a vector of all coils
    StdVector<DOMElement*> *coilSec = NULL;
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    keyVec.Push_back( "pdeList" );
    attrVec.Push_back( "" );
    valVec.Push_back( "" );
    if ( pde != "" ) {
      keyVec.Push_back( pde );
      attrVec.Push_back( "" );
      valVec.Push_back( "" );
    }
    keyVec.Push_back( "coils" );
    attrVec.Push_back( "" );
    valVec.Push_back( "" );
    coilSec = FindMatchingElements( keyVec, attrVec, valVec, rootElem_, 1 );
    if ( coilSec->GetSize() == 0 ) {
      Info->Error( "Cannot find a 'coils' section", __FILE__, __LINE__ );
    }
    else if ( coilSec->GetSize() > 1 ) {
      Info->Error( "Cannot deal with multiple 'coils' sections", __FILE__,
                   __LINE__);
    }
    DOMNodeList *coils = (*coilSec)[0]->getChildNodes();

    // Now search for coil with matching name attribute
    // and determine the element's tag, which encodes the
    // type of coil
    unsigned int foundCoils = 0;
    DOMElement *elem = NULL;
    char *auxString = NULL;

    for ( unsigned int k = 0; k < coils->getLength(); k++ ) {

      if ( coils->item(k)->getNodeType() == DOMNode::ELEMENT_NODE ) {

        elem = Node2Elem( coils->item(k) );

        if ( beVerbose_ == true ) {
          std::cerr << " coilType = " << elem->getTagName() << std::endl;
        }

        if ( AttribHasValue( elem, "name", coilName ) ) {
          auxString = X2C( elem->getTagName() );
          coilType.assign( auxString );
          FreeC( &auxString );
          foundCoils++;
        }
      }
    }

    // Check for errors
    if ( foundCoils == 0 ) {
      (*error) << "Found no coil with name '" << coilName << "'";
      Error( __FILE__, __LINE__ );
    }
    else if ( foundCoils > 1 ) {
      (*error) << "Found " << foundCoils + " coils with name '"
               << coilName << "'";
      Error( __FILE__, __LINE__ );
    }

    // Cleanup
    coilSec->Clear();
    delete coilSec;

  }


  // ========================================
  //   Query on/off status of a flag/switch
  // ========================================
  Boolean XMLParamHandler::IsSet( const std::string key,
                                  const std::string section,
                                  const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::IsSet" );

    Boolean flagStatus = FALSE;

    // Generate vectors of keywords and side-constraints
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    GenerateSearchParams( key, section, subsection, keyVec, attrVec, valVec );

    // Find all elements/values matching keyword in (restricted) tree
    StdVector<std::string> matches;
    GetList( keyVec, attrVec, valVec, matches );

    // If there is no match, return false
    if ( matches.GetSize() == 0 ) {
      flagStatus = FALSE;
    }

    // If there is a match, but it is not unique, call problem handler
    else if ( matches.GetSize() > 1 ) {
      MultipleMatchHandler( keyVec, attrVec, valVec, matches.GetSize() );
    }

    // So, there is a matching parameter. Thus, test its value
    else if ( matches[0] == "yes" ) {
      flagStatus = TRUE;
    }

    // Parameter value is not "yes"
    else {
      flagStatus = FALSE;
    }

    // we are done
    return flagStatus;

  }


  // =================================================
  //   Query whether a parameter has a certain value
  // =================================================
  Boolean XMLParamHandler::HasValue( const std::string key,
                                     const std::string value,
                                     const std::string section,
                                     const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::HasValue" );

    Boolean flagStatus = FALSE;

    // Generate vectors of keywords and side-constraints
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    GenerateSearchParams( key, section, subsection, keyVec, attrVec, valVec );

    // Find all elements/values matching keyword in (restricted) tree
    StdVector<std::string> matches;
    GetList( keyVec, attrVec, valVec, matches );

    // If there is no match, check for default
    if ( matches.GetSize() == 0 ) {
      std::string defaultValue;
      flagStatus = CheckForDefault( keyVec, attrVec, valVec, defaultValue );

      // If there is a default, then test its value.
      // If it does not match, then re-set status
      if ( flagStatus == TRUE ) {
        if ( defaultValue != value ) {
          flagStatus = FALSE;
        }
      }
    }

    // If there is a match, but it is not unique, call problem handler
    else if ( matches.GetSize() > 1 ) {
      MultipleMatchHandler( keyVec, attrVec, valVec, matches.GetSize() );
    }

    // So, there is a matching parameter. Thus, test its value
    else if ( matches[0] == value ) {
      flagStatus = TRUE;
    }

    // Parameter value is unequal to speficied value
    else {
      flagStatus = FALSE;
    }

    // we are done
    return flagStatus;
  }


  // *************************************************************************
  //   Private Auxilliary Methods: Conversion Routines
  // *************************************************************************


  // =====================================
  //   Convert a DOMNode to a DOMElement
  // =====================================
  DOMElement* XMLParamHandler::Node2Elem( DOMNode *node ) {
    ENTER_IFCN( "XMLParamHandler::Node2Elem" );
    if ( node->getNodeType() != DOMNode::ELEMENT_NODE ) {
      std::string errmsg = "Invalid conversion attempt from DOMNode to ";
      errmsg += "DOMElement!";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }
    return dynamic_cast<DOMElement *>(node); 
  }


  // ==================================
  //   Convert a DOMNode to a DOMAttr
  // ==================================
  DOMAttr* XMLParamHandler::Node2Attr( DOMNode *node ) {
    ENTER_IFCN( "XMLParamHandler::Node2Attr" );
    if ( node->getNodeType() != DOMNode::ATTRIBUTE_NODE ) {
      std::string errmsg = "Invalid conversion attempt from DOMNode to";
      errmsg += " DOMAttr!";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }
    return dynamic_cast<DOMAttr *>(node); 
  }


  // =================================
  //   Convert a DOMNode to a string
  // =================================
  void XMLParamHandler::Node2String( DOMNode *textnode, std::string &value ) {

    ENTER_IFCN( "XMLParamHandler::Node2String" );

    if ( textnode->getNodeType() != DOMNode::TEXT_NODE ) {
      (*error) << "Invalid conversion attempt from DOMNode to string!\n"
               << "Node is not of type TEXT_NODE\n Conversion of arbitrary"
               << " types to string is currently not supported!";
      Error( __FILE__, __LINE__ );
    }
    else {
      char *auxString = X2C( textnode->getNodeValue() );
      value.assign( auxString );
      FreeC( &auxString );
    }
  }


  // *************************************************************************
  //   Private Auxilliary Methods: Search Routines
  // *************************************************************************


  // ========================================
  //   Recursive search method for elements
  // ========================================
  StdVector<DOMElement *>*
  XMLParamHandler::FindMatchingElements( StdVector<std::string> &keys,
                                         StdVector<std::string> &attribs,
                                         StdVector<std::string> &aValues,
                                         DOMElement *treetop,
                                         unsigned int curdepth ) {

    ENTER_IFCN( "XMLParamHandler::FindMatchingElements" );

    DOMElement *auxElem = NULL;
    StdVector<DOMElement *> *elemvec = NULL;
    StdVector<DOMElement *> *branchTops = NULL;

    // Perform consistency checks
    if ( keys.GetSize() == 0 ) {
      Info->Error( "FindMatchingElements: Got zero keys!", __FILE__,
                   __LINE__ );
    }
    if ( curdepth == 0 || curdepth > keys.GetSize() ) {
      (*error) << "FindMatchingElements: curdepth = " << curdepth
               << ", but should be in [1," << keys.GetSize() << "]";
      Error( __FILE__, __LINE__ );
    }
    if ( attribs.GetSize() < keys.GetSize() ||
         aValues.GetSize() < keys.GetSize() ) {
      (*error) << "attribs or aValues vector too short";
      Error( __FILE__, __LINE__ );
    }

    // If desired be verbose
    if ( beVerbose_ == true ) {
      std::cerr << " FindMatchingElements: curdepth = " << curdepth
                << "\n                       key[" << curdepth << "] = "
                << keys[curdepth-1]
                << "\n                       attrib[" << curdepth << "] = "
                << attribs[curdepth-1]
                << "\n                       aValue[" << curdepth << "] = "
                << aValues[curdepth-1]
                << std::endl;
    }

    // Generate list of all elements in subtree starting at treetop
    // that match the curdepth'th keyword
    XMLCh *tmpString = S2X( keys[curdepth-1] );
    DOMNodeList *list = treetop->getElementsByTagName( tmpString );
    FreeX( &tmpString );

    // Test, if there are side constraints for the attributes of the elements
    // on this search depth, and, if there are, remove from the list of
    // matching elements those that do not meet those side constraints
    branchTops = new StdVector<DOMElement *>;
    if ( attribs[curdepth-1] == "" ) {
      for ( unsigned int i = 0; i < list->getLength(); i++ ) {
        branchTops->Push_back( Node2Elem( list->item(i) ) );
      }
    }
    else {
      if ( beVerbose_ == true ) {
        std::cerr << " -> Testing " << list->getLength() << " element(s) that"
                  << " are potential matches" << std::endl;
      }
      for ( unsigned int i = 0; i < list->getLength(); i++ ) {
        auxElem = Node2Elem( list->item(i) );
        if ( AttribHasValue( auxElem, attribs[curdepth-1],
                             aValues[curdepth-1], false ) ) {
          branchTops->Push_back( auxElem );
        }
      }
    }

    // If we are on the lowest level allowed, we can simply return the vector
    // of matching elements
    if ( curdepth == keys.GetSize() ) {
      elemvec = branchTops;
      if ( beVerbose_ == true ) {
        std::cerr << " Found " << elemvec->GetSize() << " matches!\n";
        std::cerr << " Reached bottom of recursion!\n";
      }
    }

    // There were no new branchTops so return immediately
    else if ( branchTops->GetSize() == 0 ) {
      if ( beVerbose_ == true ) {
        std::cerr << " No new sub-tress -> ascending again\n";
      }
      elemvec = branchTops;
    }

    // This is not the lowest level, so we descend and look for elements
    // in all subtrees spanned by the matches in our element vector and
    // append the results to the element vector
    else {

      // Report descend
      if ( beVerbose_ == true ) {
        std::cerr << " Descending: Got " << branchTops->GetSize()
                  << " new subtrees\n";
      }

      // Generate results vector for this level
      elemvec = new StdVector<DOMElement *>;

      // auxilliary vector for intermediate results
      StdVector<DOMElement *> *tmpvec = NULL;

      // loop over all allowed subtrees
      for ( unsigned int i = 0; i < branchTops->GetSize(); i++ ) {

        // get results for each subtree
        tmpvec = FindMatchingElements( keys, attribs, aValues,
                                       (*branchTops)[i], curdepth+1 );

        // Report ascend
        if ( beVerbose_ == true ) {
          std::cerr << " Back on level " << curdepth << '\n';
        }

        // append results to our element vector
        for ( unsigned int k = 0; k < tmpvec->GetSize(); k++ ) {
          elemvec->Push_back( (*tmpvec)[k] );
        }

        // Clean up
        tmpvec->Clear();
        delete tmpvec;
      }

      // At end of ascending we delete intermediate stuff that is no
      // longer necessary to avoid memory leaks and _yes_, it _is_
      // intended that this delete is not deep!
      delete branchTops;
    }

    // Return element vector
    return elemvec;
  }


  // ==========================================
  //   Recursive search method for attributes
  // ==========================================
  StdVector<DOMAttr*>*
  XMLParamHandler::FindMatchingAttributes( std::string attr_key,
                                           StdVector<std::string> &keys,
                                           StdVector<std::string> &attribs,
                                           StdVector<std::string> &aValues,
                                           DOMElement *treeTop ) {

    ENTER_IFCN( "XMLParamHandler::FindMatchingAttributes" );

    // Perform consistency checks
    if ( keys.GetSize() == 0 ) {
      Info->Error( "FindMatchingAttributes: Got zero keys!", __FILE__,
                   __LINE__ );
    }

    // Generate return vector
    StdVector<DOMAttr*> *attrvec = new StdVector<DOMAttr*>;

    // Find all matching elements
    StdVector<DOMElement *>* elements =
      FindMatchingElements( keys, attribs, aValues, treeTop, 1 );

    // Be verbose, if demanded
    if ( beVerbose_ == true ) {
      std::cerr << " FindMatchingAttributes: Found " << elements->GetSize()
                << " matches" << std::endl;
    }

    // Loop over all elements and examine their attributes
    DOMNamedNodeMap *attributes = NULL;
    DOMNode* aux_node = NULL;
    DOMAttr* match_attr = NULL;
    std::string dummy;
    XMLCh *tmpString = NULL;

    for ( unsigned int i = 0; i < elements->GetSize(); i++ ) {

      // Get attribute of elements
      attributes = (*elements)[i]->getAttributes();

      // Select attributes matching name and convert them
      tmpString = S2X( attr_key );
      aux_node = attributes->getNamedItem( tmpString );
      FreeX( &tmpString );
      if ( aux_node != NULL ) {
        match_attr = Node2Attr( aux_node );
      }
      else {
        match_attr = NULL;
      }
      
      // If there is a match, append attribute to result vector
      if ( match_attr != NULL ) {
        attrvec->Push_back( match_attr );
      }
    }

    // Free space allocated by ourselves in this routine
    if ( elements != NULL ) {
      elements->Clear();
      delete elements;
    }

    // Finished
    return attrvec;

  }


  // ===============================================
  //   Return a list of strings matching a keyword
  // ===============================================
  void XMLParamHandler::FindAllMatches( const StdVector<std::string>& key,
                                        const StdVector<std::string>& attrib,
                                        const StdVector<std::string>& aValue,
                                        StdVector<std::string> &list,
                                        DOMElement *treeTop ) {

    ENTER_IFCN( "XMLParamHandler::FindAllMatches" );

    // *************************************
    //   Part 1: Check of input parameters
    // *************************************

    // Perform consistency checks
    if ( key.GetSize() == 0 ) {
      Info->Error( "FindAllMatches: Got zero keys!", __FILE__, __LINE__ );
    }
    if ( attrib.GetSize() != key.GetSize()-1 ||
         aValue.GetSize() != key.GetSize()-1 ) {
      std::cerr << "key vector: " << key.GetSize() << '\n'
                << "attribute vector: " << attrib.GetSize() << '\n'
                << "value vector: " << aValue.GetSize() << '\n';
      Info->Error( "Improper length of attribs or aValues vector", __FILE__,
                   __LINE__ );
    }

    // Report, what we are looking for
    if ( beVerbose_ == true ) {
      fprintf( stderr, "\n FindAllMatches: Starting new search\n" );
      PrintSearchParams( key, attrib, aValue, stderr );
    }

    // Check if vector is empty. If not issue a warning
    // and erase its entries, if this is desired
    ClearVector( list );


    // **************************
    //   Part 2: Prepare search
    // **************************
    StdVector<std::string> myKeys       = StdVector<std::string>( key    );
    StdVector<std::string> attribNames  = StdVector<std::string>( attrib );
    StdVector<std::string> attribValues = StdVector<std::string>( aValue );
    attribNames.Push_back( "" );
    attribValues.Push_back( "" );

    bool elem_match = false;
    bool attr_match = false;
    StdVector<DOMElement*>* elem_matches = NULL;
    StdVector<DOMAttr*   >* attr_matches = NULL;
    StdVector<DOMElement*> attr_elements;
    StdVector<std::string> keys;
    std::string value;


    // **************************************************
    //   Part 3: Assume main key word is an element tag
    // **************************************************

    // Find matching elements
    elem_matches = FindMatchingElements( myKeys, attribNames, attribValues,
                                         treeTop, 1 );

    // Check if there was a match
    if ( elem_matches->GetSize() > 0 ) {
      elem_match = true;
    }

    // Be verbose, if demanded
    if ( beVerbose_ == true ) {
      std::cerr << " FindAllMatches: Found " << elem_matches->GetSize()
                << " matching elements" << std::endl;
    }

    // *****************************************************
    //   Part 4: Assume main key word is an attribute name
    // *****************************************************

    // Test whether a section was specified. If not, no attribute
    // search is possible
    unsigned int numKeys = myKeys.GetSize();
    if ( numKeys >= 2 ) {

      // Last key is the name of the attribute
      std::string attr_key = myKeys[numKeys-1];

      // Shorten vectors, since only the first to last but one entry
      // represents the name of an element
      myKeys.Erase( numKeys-1 );
      attribNames.Erase( numKeys-1 );
      attribValues.Erase( numKeys-1 );

      // Find matching attributes
      attr_matches = FindMatchingAttributes( attr_key, myKeys, attribNames,
                                             attribValues, treeTop );

      // Be verbose, if demanded
      if ( beVerbose_ == true ) {
        std::cerr << " FindAllMatches: Found " << attr_matches->GetSize()
                  << " matching attributes" << std::endl;
      }

      // Check if there was a match, if not call problem handler
      if ( attr_matches->GetSize() > 0 ) {
        attr_match = true;
      }
    }

    // **********************************
    //   Part 5: Process search results
    // **********************************

    char *auxString = NULL;

    // See, if there were element and attribute matches
    if ( elem_match == true && attr_match == true ) {
      (*error) << "Keyword '" << key[key.GetSize()-1]
               << "' matches element(s) and attribute(s)!";
      Error( __FILE__, __LINE__ );
    }

    // Convert element values to strings and append to list
    else if ( elem_match == true ) {
      for ( unsigned int i = 0; i < elem_matches->GetSize(); i++ ) {
        auxString = GetElementValue( (*elem_matches)[i] );
        value.assign( auxString );
        list.Push_back( value );
        FreeC( &auxString );
      }
    }

    // Convert element values to strings and append to list
    else if ( attr_match == true ) {
      for ( unsigned int i = 0; i < attr_matches->GetSize(); i++ ) {
        auxString = X2C( (*attr_matches)[i]->getValue() );
        value.assign( auxString );
        list.Push_back( value );
        FreeC( &auxString );
      }
    }

    // *******************
    //   Part 6: Report
    // *******************
    if ( beVerbose_ == true ) {
      std::cerr << "\n FindAllMatches found the following matches:\n";
      for ( unsigned int i = 0; i < list.GetSize(); i++ ) {
        std::cerr << "match[" << i << "] = '" << list[i] << "'\n\n";
      }
    }

    // *******************
    //   Part 7: Cleanup
    // *******************
    if ( elem_matches != NULL ) {
      elem_matches->Clear();
      delete elem_matches;
    }
    if ( attr_matches != NULL ) {
      attr_matches->Clear();
      delete attr_matches;
    }
  }


  // *************************************************************************
  //   Private Auxilliary Methods: Treatment of errors and defaults
  // *************************************************************************


  // ==================================
  //   Treat case of multiple matches
  // ==================================
  void XMLParamHandler::
  MultipleMatchHandler( const StdVector<std::string> &keyVec,
                        const StdVector<std::string> &attrVec,
                        const StdVector<std::string> &valVec,
                        const unsigned int nmatches ) {

    ENTER_IFCN( "XMLParamHandler::MultipleMatchHandler" );

    std::string errmsg;

    // Make sure that we have a non-unique match
    if ( nmatches > 1 ) {
      fprintf( stderr, "\n\n XMLParamHandler:\n Found %d matches ",
               nmatches );
      fprintf( stderr, " while searching for:\n\n" );
      PrintSearchParams( keyVec, attrVec, valVec, stderr );

      std::string errmsg;
      errmsg += "No unique match found!";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // This should not happen!
    else {
      errmsg  = "Internal Error: MultipleMatchHandler called with ";
      errmsg += "nmatches = ";
      errmsg += Info->GenStr(nmatches);
      errmsg += "!\n This should not have happened!\n";
    }

    // And now terminate
    Info->Error( errmsg, __FILE__, __LINE__ );

  }


  // ==========================
  //   Treat case of no match
  // ==========================
  void XMLParamHandler::NoMatchHandler( const StdVector<std::string> &keyVec,
                                        const StdVector<std::string> &attrVec,
                                        const StdVector<std::string> &valVec,
                                        std::string &defaultValue ) {

    ENTER_IFCN( "XMLParamHandler::NoMatchHandler" );

    // Test, whether a default value is specified for the parameter
    Boolean defaultExists = CheckForDefault( keyVec, attrVec, valVec,
                                             defaultValue );

    // If no default could be found, cry out!
    if( defaultExists == FALSE ) {
      NoMatchErrorReporter( keyVec, attrVec, valVec );
    }
  }


  // =================================
  //   Report that there is no match
  // =================================
  void XMLParamHandler::
  NoMatchErrorReporter( const StdVector<std::string> &keyVec,
                        const StdVector<std::string> &attrVec,
                        const StdVector<std::string> &valVec ) {

    ENTER_IFCN( "XMLParamHandler::NoMatchErrorReporter" );

    fprintf( stderr, "\n\n XMLParamHandler:\n No match and no default" );
    fprintf( stderr, " found while searching for:\n\n" );
    PrintSearchParams( keyVec, attrVec, valVec, stderr );

    std::string errmsg;
    errmsg += "No match and no default found!";
    Info->Error( errmsg, __FILE__, __LINE__ );
  }


  // =================================
  //   Check for default parameter
  // =================================
  Boolean
  XMLParamHandler::CheckForDefault( const StdVector<std::string> &keyVec,
                                    const StdVector<std::string> &attrVec,
                                    const StdVector<std::string> &valVec,
                                    std::string &defaultValue ) {

    ENTER_IFCN( "XMLParamHandler::CheckForDefault" );

    Boolean defaultFound = FALSE;
    StdVector<std::string> newValVec( valVec );

    // Check if string is empty. If not issue a warning
    // and erase it, if this is desired
    if ( defaultValue.size() != 0 ) {
      if ( beVerbose_ == true ) {
        std::string msg = "Warning input string was not empty!\n";
        msg += "Contents have been erased!";
        Info->Warning( msg );
      }
      defaultValue.clear();
    }

    // NOTE: There is a conceptual problem with a restricted search and a
    //       default parameter tree, e.g. we can look for a property of a
    //       coil and select that coil by restricting the name attribute,
    //       however in the default XML file we cannot now that name.
    //       This is a brute force attempt to fix this for the above coil
    //       example.
    for ( unsigned int i = 0; i < keyVec.GetSize(); i++ ) {
      if ( keyVec[i] == "coils" && attrVec[i+1] == "name" ) {
        newValVec[i+1] = "dummyRegion";
        if ( beVerbose_ == true ) {
          Info->Warning( "Adapted name of coil to dummyRegion" );
        }
        break;
      }
    }

    // Search for matching elements/attributes in default tree
    StdVector<std::string> matches;
    FindAllMatches( keyVec, attrVec, newValVec, matches, rootElemDefaults_ );

    // If there was no unique match, call problem handler
    if ( matches.GetSize() > 1 ) {

      // Test if matches are different
      bool valsAgree = true;
      for ( UInt k = 1; k < matches.GetSize(); k++ ) {
        if ( matches[k] != matches[0] ) {
          valsAgree = false;
        }
      }

      // values do not agree
      if ( valsAgree == false ) {
        MultipleMatchHandler( keyVec, attrVec, valVec, matches.GetSize() );
      }

      // values agree
      else {
        defaultFound = TRUE;
        defaultValue = matches[0];
      }
    }

    // Check, if a default was found
    if ( matches.GetSize() == 1 ) {
      defaultFound = TRUE;
      defaultValue = matches[0];
    }

    // If no default was found, test, whether this is one of the special
    // cases not implemented so far
    else {

      if( keyVec[keyVec.GetSize()-1] == "preStressVal" ) {
        defaultFound = TRUE;
        defaultValue = "0";
      }

      if( keyVec[keyVec.GetSize()-1] == "effMass" ) {
        defaultFound = TRUE;
        defaultValue = "no";
      }

    }

    // Tell what we found
    if ( beVerbose_ == true && defaultFound == TRUE ) {
      std::string msg = "CheckForDefault: Default for parameter '" +
        keyVec[keyVec.GetSize()-1];
      msg += "' is '" + defaultValue + "'";
      Info->Warning( msg );
    }

    return defaultFound;

  }


  // *************************************************************************
  //   Private Auxilliary Methods: Miscellaneous
  // *************************************************************************


  // ===========================
  //   Get value of an element
  // ===========================
  char* XMLParamHandler::GetElementValue( DOMElement *elem ) {

    ENTER_IFCN( "XMLParamHandler::GetElementValue" );

    // Obtain child nodes of this element
    DOMNodeList *children = elem->getChildNodes();

    // Complain, if there is not exactly one child
    if ( children->getLength() != 1 ) {
      if ( children->getLength() == 0 ) {
        (*error) << "GetElementValue: Encountered element without child!";
      }
      else {
        (*error) << "GetElementValue: Encountered element with multiple "
                 << "children!";
      }
      (*error) << "\n       GetElementValue: Element tag is '"
               << elem->getNodeName() << "'";
      Error( __FILE__, __LINE__ );
    }

    // Make sure that child is a text node
    DOMNode *child = children->item(0);
    if ( child->getNodeType() != DOMNode::TEXT_NODE ) {
      Info->Error( "GetElementValue: Child of element is no TEXT_NODE!",
                   __FILE__, __LINE__ );
    }
    
    // Convert value of element to character array
    return X2C( child->getNodeValue() );
  }


  // =======================================
  //   Get value of an element's attribute
  // =======================================
  Boolean XMLParamHandler::GetElementAttribute( DOMElement* element,
                                                const std::string keyword,
                                                std::string &attrVal ) {

    ENTER_IFCN( "XMLParamHandler::GetElementAttribute" );

    bool hasAttr;
    DOMNamedNodeMap *attributes = NULL;
    DOMNode* aux_node = NULL;
    DOMAttr* aux_attr = NULL;
    char *auxString = NULL;

    // Get attribute of elements
    attributes = element->getAttributes();

    // Select attributes matching name and convert them
    XMLCh *tmpString = S2X( keyword );
    aux_node = attributes->getNamedItem( tmpString );
    FreeX( &tmpString );

    // Check, if element has specified attribute
    if ( aux_node != NULL ) {
      hasAttr = true;
    }
    else {
      hasAttr = false;
    }

    // If element has specified attribute obtain its value
    if ( hasAttr ) {
      aux_attr = Node2Attr( aux_node );
      auxString = X2C( aux_attr->getValue() );
      attrVal.assign( auxString );
      FreeC( &auxString );
    }
    else {
      attrVal = "";
    }

    // Finished
    if ( hasAttr ) {
      return TRUE;
    }
    return FALSE;
  }


  // ============================
  //   Parse XML parameter file
  // ============================
  DOMElement* XMLParamHandler::ParseFile( XercesDOMParser **parser,
                                          const char *xmlFile ) {

    // Root element of DOMTree built by parser
    DOMElement *rootElem = NULL;

    // Create the parser, force validation using full schema checking with
    // namespaces and tell it to drop irrelevant white space
    (*parser) = new XercesDOMParser();
    (*parser)->setValidationScheme(XercesDOMParser::Val_Always);
    (*parser)->setDoNamespaces(true);
    (*parser)->setDoSchema(true);
    (*parser)->setValidationSchemaFullChecking(true);
    (*parser)->setIncludeIgnorableWhitespace(false);

    // We may separate the schema file from the instance file
    cfsSchema_ = "http://www.cfs++.org ";
    cfsSchema_ += commandLine->GetSchemaPath();
    cfsSchema_ += "/CFS.xsd";
    (*parser)->setExternalSchemaLocation( cfsSchema_.c_str() );
    
    // Have not yet understood what an entity reference node is,
    // but it seems we do not need them
    (*parser)->setCreateEntityReferenceNodes(false);

    // Attach our own error handler to the parser
    XMLParserErrorHandler *errHandler = new XMLParserErrorHandler();
    (*parser)->setErrorHandler(errHandler);

    // Parse and validate the XML file. This will generate the DOM tree.
    // Catch all exceptions that the parser could not pass to our error
    // handler.
    try {
      (*parser)->parse( xmlFile );
    }
    catch( const XMLException &event ) {
      (*error) << "The following XML error was encountered during "
               << "parsing:\n" << event.getMessage();
      Error( __FILE__, __LINE__ );
    }
    catch( const DOMException &event ) {
      const unsigned int maxChars = 2047;
      XMLCh errText[maxChars + 1];

      (*error) << "DOM Error during parsing! DOMException code is:\n"
               << event.code << '\n';

      if( DOMImplementation::loadDOMExceptionMsg(event.code, errText,
                                                 maxChars) ) {
        (*error) << "Message is: " << errText;
      }

      Error( __FILE__, __LINE__ );
    }
    catch(...) {
      (*error) << "An unknown error occurred during parsing!\n"
               << "All I can say is that it was neither an XMLException "
               << "nor a DOMException.";
      Error( __FILE__, __LINE__ );
    }

    // Obtain and validate root element of document tree
    DOMDocument *doc = (*parser)->getDocument();
    DOMNodeList *children = doc->getChildNodes();
    if ( children->getLength() != 1 ) {
      (*error) << "Document root has more than one child!\n"
               << "We are in real trouble here!";
      Error( __FILE__, __LINE__ );
    }
    if ( children->item(0)->getNodeType() != DOMNode::ELEMENT_NODE ) {
      (*error) << "Root node is not of type DOMNode::ELEMENT_NODE!\n"
               << "We are in real trouble here!";
      Error( __FILE__, __LINE__ );
    }
    rootElem = (DOMElement *)(children->item(0)); 

    // We are finished
    return rootElem;

  }


  // =====================================================
  //   Check if element has attribute with certain value
  // =====================================================
  bool XMLParamHandler::AttribHasValue( DOMElement* elem,
                                        const std::string attribute,
                                        const std::string value,
                                        bool failIfNoAttrib ) {

    ENTER_IFCN( "XMLParamHandler::AttribHasValue" );

    bool retVal = false;

    // Get element's attributes
    DOMNamedNodeMap *attributes = NULL;
    attributes = elem->getAttributes();

    // Select attribute matching name
    XMLCh *tmpString = S2X( attribute );
    DOMNode* attrib = attributes->getNamedItem( tmpString );
    FreeX( &tmpString );

    StdVector<std::string> attribs;

    // Test, if element has an attribute with specified name
    if ( attrib == NULL ) {
      if ( failIfNoAttrib == true ) {
        (*error) << "AttribHasValue: Element '" << elem->getTagName()
                 << "' has no attribute '" << attribute << "'";
        Error( __FILE__, __LINE__ );
      }
      else {
        retVal = false;
      }
    }

    // Test, if attribute's value matches specification
    else {

      // Report for debugging
      if ( beVerbose_ == true ) {
        std::cerr << " Checking attribute '" << attribute
                  << "' for value '" << value << "'" << std::endl;
      }

      // Convert attribute's value to string representation
      std::string *foundValue = X2S( Node2Attr(attrib)->getValue() );

      // Treat case of 'tag' attribute in special way
      if ( attribute == "tag" ) {
        retVal = MatchesTag( *foundValue, value );
      }

      // Standard case
      else {
        retVal = *foundValue == value ? true : false;
      }

      FreeS( &foundValue );
    }

    return retVal;
  }
  

  // ===============================
  //   Check if refTag matches tag
  // ===============================
  bool XMLParamHandler::MatchesTag( const std::string tag,
                                    const std::string refTag ) {

    ENTER_IFCN( "XMLParamHandler::AttribHasValue" );

    bool match = false;

    // Case 1: tag matches always
    if ( tag == "anyTag" ) {
      match = true;
    }

    // Case 2: one-to-one match
    else if ( tag == refTag ) {
      match = true;
    }

    // Case 3: refTag contained as sub-tag
    else {
      StdVector<std::string> subTags;
      SplitStringList( tag, subTags, ',' );
      for ( unsigned int i = 0; i < subTags.GetSize(); i++ ) {
        if ( subTags[i] == refTag ) {
          match = true;
          break;
        }
      }
    }

    // done
    return match;

  }


  // ===============================================
  //   Print search parameters to an output stream
  // ===============================================
  void
  XMLParamHandler::PrintSearchParams( const StdVector<std::string> &key,
                                      const StdVector<std::string> &attr,
                                      const StdVector<std::string> &value,
                                      FILE *myStream ) {

    ENTER_IFCN( "XMLParamHandler::PrintSearchParams" );

    fprintf( myStream, "  key\t\tattribute\tvalue\n" );
    for ( unsigned int i = 0; i < key.GetSize() - 1; i++ ) {
      fprintf( myStream, "  '%s'\t'%s'\t'%s'\n", key[i].c_str(),
               attr[i].c_str(), value[i].c_str() );
    }
    fprintf( myStream, "  '%s'\t''\t''\n", key[key.GetSize()-1].c_str() );

  }



  // =====================================================
  //   Generate vectors of keywords and side constraints
  // =====================================================
  void 
  XMLParamHandler::GenerateSearchParams( const std::string key,
                                         const std::string section,
                                         const std::string subsection,
                                         StdVector<std::string> &keyVec,
                                         StdVector<std::string> &attrVec,
                                         StdVector<std::string> &valVec ) {

    ENTER_FCN( "XMLParamHandler::GenerateSearchParams" );

    // Check for section
    if ( section != "" ) {
      keyVec.Push_back( section );
      attrVec.Push_back( "" );
      valVec.Push_back( "" );
    }

    // Check for subsection
    if ( subsection != "" ) {
      keyVec.Push_back( subsection );
      attrVec.Push_back( "" );
      valVec.Push_back( "" );
    }

    // There should be a key
    keyVec.Push_back( key );
  }


  // =====================================
  //   Print DOM tree to standard output
  // =====================================
  void XMLParamHandler::PrintTree( xercesc::DOMDocument *doc ) {

    ENTER_FCN( "XMLParamHandler::PrintTree" );

    try {

      // Obtain a DOMImplementationLS implementation from the registry and
      // use this factory to generate an instance of DOMWriter
      XMLCh *tmpString = C2X( "LS" );
      DOMImplementation *impl =
        DOMImplementationRegistry::getDOMImplementation( tmpString );
      FreeX( &tmpString );
      DOMWriter *theSerializer =
        ((DOMImplementationLS*)impl)->createDOMWriter();

      // Plug in a format target to receive the resultant
      // XML stream from the serializer.
      //
      // StdOutFormatTarget prints the resultant XML stream
      // to stdout once it receives anything from the serializer.
      XMLFormatTarget *myFormTarget = new StdOutFormatTarget();

      // Set some features
      DOMWriterSetFeature( theSerializer, XMLUni::fgDOMWRTFormatPrettyPrint,
                           true, true );
      DOMWriterSetFeature( theSerializer,
                           XMLUni::fgDOMWRTDiscardDefaultContent,
                           false, true );

      // Report what comes next
      std::cout << "\n"
                << " --------------------------------------------------------"
                << "----------------------\n"
                << "   The following is a serialisation of the DOM tree "
                << "generated by the parser.\n   It includes all values "
                << "for optional attributes that were not specified in\n"
                << "   the XML file and were added by the parser based on the"
                << " Schema definitions\n   used for validation.\n\n"
                << "   Note:\n\n"
                << "   The W3C Schema specifications demand that only "
                << " default values for un-\n   supplied attributes"
                << " are added by the parser. Optional elements are not\n"
                << "   added to the tree and default values are only"
                << " supplied for empty but\n   not for missing elements.\n"
                << " --------------------------------------------------------"
                << "----------------------\n"
                << std::endl;

      // Perform the serialization through DOMWriter::writeNode();
      theSerializer->writeNode( myFormTarget, *doc );

      // Clean up
      delete theSerializer;

    }
    catch (XMLException& e) {
      std::cerr << "Some weird shit happened!\n" << std::endl;
    }
  }


  // ================================
  //   Set feature of the DOMWriter
  // ================================
  bool XMLParamHandler::DOMWriterSetFeature( DOMWriter *serialiser,
                                             const XMLCh *const feature,
                                             bool featureVal,
                                             bool shouldHave ) {

    ENTER_FCN( "XMLParamHandler::DOMWriterSetFeature" );

    if ( featureVal == true ) {
      std::cout << " Enabling feature '";
    }
    else {
      std::cout << " Disabling feature '";
    }
    std::cout << feature << "'." << std::endl;

    // Check if serialiser supports setting the feature
    // to the specified value
    bool supportsFeature = serialiser->canSetFeature( feature, featureVal );
    if ( supportsFeature == false ) {
      if ( shouldHave == true ) {
        std::cout << "\n ERROR: The serialiser does not support enabling "
                  << "feature '" << feature << "'\n\n";
        exit(-1);
      }
      else {
        std::cout << "\n WARNING: The serialiser does not support enabling "
                  << "feature '" << feature << "'\n\n";
      }
    }

    // The serialiser supports it, so enable/disable the feature
    else {
      serialiser->setFeature( feature, featureVal );
    }

    // Return, whether the serialiser supported setting
    // the feature
    return supportsFeature;
  }

}

#endif
