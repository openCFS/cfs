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

#include "General/environment.hh"
#include "DataInOut/WriteInfo.hh"
#include "BaseParamHandler.hh"
#include "XMLParamHandler.hh"
#include "XMLParserErrorHandler.hh"

// For internal debugging
#define LOGTOINFO false

namespace CoupledField {


  // **************************************************************************
  //   Constructors and Destructors
  // **************************************************************************


  // =========================================
  //   Constructor for XMLParamHandler class
  // =========================================
  XMLParamHandler::XMLParamHandler( char *fname ): parser_(NULL) {

    ENTER_FCN( "XMLParamHandler::XMLParamHandler" );

    // String for assembling error messages
    std::string errmsg;

    // Initialise the XML4C2 system
    try {
      XMLPlatformUtils::Initialize();
    }
    catch( const XMLException &event ) {
      errmsg  = "The following error occured during initialisation of ";
      errmsg += "xerces-c!\n";
      errmsg += X2S(event.getMessage());
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // Generate parser and parse XML parameter file
    rootElem_ = ParseFile( &parser_, fname );

    // Generate parser and parse XML defaults file
    std::string cfsDefaults = CVSEXTERNAL;
    cfsDefaults += "/CFS++XML/Defaults/CFS++Defaults.xml";
    rootElemDefaults_ = ParseFile( &parserDefaults_, cfsDefaults.c_str() );

    // Toggle verbosity
#ifdef DEBUG_XMLPARAMHANDLER
    beVerbose_ = true;
#else
    beVerbose_ = false;
#endif

  }


  // ========================================
  //   Destructor for XMLParamHandler class
  // ========================================
  XMLParamHandler::~XMLParamHandler() {
    ENTER_FCN( "XMLParamHandler::~XMLParamHandler" );
    if ( parser_ != NULL ) {
      delete parser_;
    }
    if ( parserDefaults_ != NULL ) {
      delete parserDefaults_;
    }
  }


  // =================================================
  //   Default constructor for XMLParamHandler class
  // =================================================
  XMLParamHandler::XMLParamHandler() {}


  // **************************************************************************
  //   Public Methods
  // **************************************************************************


  // ================================================================
  //   Return as string the value of a parameter matching a keyword
  // ================================================================
  void XMLParamHandler::Get( const std::string key, std::string &value,
			     const std::string section,
			     const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::Get" );

    // Find all elements/values matching keyword in (restricted) tree
    std::vector<std::string> matches;
    GetList( key, matches, section, subsection );

    // If there was no unique match, call problem handler
    if ( matches.size() > 1 ) {
      MultipleMatchHandler( key, section, subsection, matches.size() );
    }

    // If there was no match at all, call problem handler
    else if ( matches.size() == 0 ) {
      NoMatchHandler( value, key, section, subsection );
    }
    
    // There was a unique match, so convert detected value
    else {
      value = matches[0];
    }

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

    // Find all elements/values matching keyword in (restricted) tree
    std::vector<std::string> matches;
    GetList( key, matches, section, subsection );

    // If there was no unique match, call problem handler
    if ( matches.size() > 1 ) {
      MultipleMatchHandler( key, section, subsection, matches.size() );
    }

    // If there was no match at all, call problem handler
    else if ( matches.size() == 0 ) {
      NoMatchHandler( value, key, section, subsection );
    }
    
    // There was a unique match, so convert detected value
    else {
      const char *val = matches[0].c_str();
      value = atoi(val);
    }

    // Tell what we found
    if ( beVerbose_ == true ) {
	std::string msg = "Get: Value for parameter '" + key;
	msg += "' is '" + Info->GenStr(value) + "'";
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

    // Find all elements/values matching keyword in (restricted) tree
    std::vector<std::string> matches;
    GetList( key, matches, section, subsection );

    // If there was no unique match, call problem handler
    if ( matches.size() > 1 ) {
      MultipleMatchHandler( key, section, subsection, matches.size() );
    }

    // If there was no match at all, call problem handler
    else if ( matches.size() == 0 ) {
      NoMatchHandler( value, key, section, subsection );
    }
    
    // There was a unique match, so convert detected value
    else {
      const char *val = matches[0].c_str();
      value = atof(val);
    }

    // Tell what we found
    if ( beVerbose_ == true ) {
	std::string msg = "Get: Value for parameter '" + key;
	msg += "' is '" + Info->GenStr(value) + "'";
	Info->Warning( msg );
    }

  }


  // ====================================================
  //   Return list of strings values matching a keyword
  // ====================================================
  void XMLParamHandler::GetList( const std::string key,
				 std::vector<std::string> &list,
				 const std::string section,
				 const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::GetList" );

    // This is easy, we just have to pass everything through to the
    // central auxilliary search method
    FindAllMatches( key, list, section, subsection, rootElem_ );

  }


  // ===================================================
  //   Return list of Double values matching a keyword
  // ===================================================
  void XMLParamHandler::GetList( const std::string key,
				 std::vector<Double> &list,
				 const std::string section,
				 const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::GetList" );

    // Check if vector is empty. If not issue a warning
    // and erase its entries, if this is desired
    if ( list.empty() != true ) {
      if ( beVerbose_ == true ) {
	std::string errmsg  = "Warning input vector was not empty!\n";
	errmsg += "Contents have been erased!";
	Info->Warning( errmsg );
      }
      list.clear();
    }

    // First determine all matches as strings
    std::vector<std::string> matches;
    FindAllMatches( key, matches, section, subsection, rootElem_ );

    // Now perform conversion
    for ( unsigned int i = 0; i < matches.size(); i++ ) {
      list.push_back( String2Double( matches[i] ) );
    }

  }


  // =====================================
  //   Return a list of the defined PDEs
  // =====================================
  void XMLParamHandler::GetPDEList( std::vector<std::string> &list ) {

    ENTER_FCN( "XMLParamHandler::GetPDEList" );

    // string for assembling error messages
    std::string errmsg;

    // Check if vector is empty. If not issue a warning
    // and erase its entries, if this is desired
    if ( list.empty() != true ) {
      if ( beVerbose_ == true ) {
	errmsg  = "Warning input vector was not empty!\n";
	errmsg += "Contents have been erased!";
	Info->Warning( errmsg );
      }
      list.clear();
    }

    // Find PDE section
    DOMNodeList *pdesec = rootElem_->getElementsByTagName( C2X("pdeList") );

    // Check that there is only one such section
    if ( pdesec->getLength() != 1 ) {
      errmsg  = "Got " + Info->GenStr( pdesec->getLength() );
      errmsg += " pdeList elements in parameter file!";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // The names of the PDEs are the tags of the child elements of the
    // PDE_list element, so get hold of them and make sure, there is
    // at least one child/PDE
    DOMNodeList *pdelist = pdesec->item(0)->getChildNodes();
    if ( pdelist->getLength() == 0 ) {
      errmsg = "Cannot find a single PDE in parameter file!";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // Now get hold of tags, convert them to strings and assemble vector
    std::string pdename;
    for ( unsigned int i = 0; i < pdelist->getLength(); i++ ) {

      // Only treat element children and not comments!
      if ( pdelist->item(i)->getNodeType() == DOMNode::ELEMENT_NODE ) {
	pdename.assign( X2C( Node2Elem( pdelist->item(i) )->getNodeName() ) );
	list.push_back( pdename );
      }
    }
  }


  // ===============================================
  //   Obtain list of attribute values for matches
  // ===============================================
  void XMLParamHandler::GetValsForHits( const std::string attribute2,
					std::vector<std::string> &vals,
					const std::string attribute1,
					const std::string keyword,
					const std::string section,
					const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::GetValsForHits" );

    // Check if vector is empty. If not issue a warning
    // and erase its entries, if this is desired
    if ( vals.empty() != true ) {
      if ( beVerbose_ == true ) {
	std::string errmsg  = "Warning input vector was not empty!\n";
	errmsg += "Contents have been erased!";
	Info->Warning( errmsg );
      }
      vals.clear();
    }

    //  Find all elements with attributes of type attribute1
    std::vector<DOMElement*> elems;
    std::vector<std::string> attrVals1;
    FindAllMatches( attribute1, attrVals1, section, subsection, rootElem_,
		    &elems );

    // Report
    if( beVerbose_ ) {
      std::string msg = "Found " + Info->GenStr( elems.size() )
	+ " elements ";
      msg += "with attribute '" + attribute1 + "' in subsection '";
      msg += subsection + "' of section '" + section + "'";
      Info->Warning( msg );
    }
    
    // Find those elements whose attribute has the correct value
    std::vector<DOMElement*> elemvec;
    for ( unsigned int k = 0; k < elems.size(); k++ ) {
      if ( attrVals1[k] == keyword ) {
	elemvec.push_back( elems[k] );
      }
    }

    // If there is no match return an empty list
    if ( elemvec.size() == 0 ) {
      return;
    }

    // For each remaining element determine value of attribute2
    std::string attrVal;
    for ( unsigned int k = 0; k < elemvec.size(); k++ ) {
      if ( GetElementAttribute( elemvec[k], attribute2, attrVal ) == TRUE ) {
	vals.push_back( attrVal );
      }
    }
  }


  // ========================================
  //   Query on/off status of a flag/switch
  // ========================================
  Boolean XMLParamHandler::IsSet( const std::string key,
				  const std::string section,
				  const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::IsSet" );

    Boolean flagStatus = FALSE;

    // Find all elements/values matching keyword in (restricted) tree
    std::vector<std::string> matches;
    GetList( key, matches, section, subsection );

    // If there is no match, return false
    if ( matches.size() == 0 ) {
      flagStatus = FALSE;
    }

    // If there is a match, but it is not unique, call problem handler
    else if ( matches.size() > 1 ) {
      MultipleMatchHandler( key, section, subsection, matches.size() );
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

    // Find all elements/values matching keyword in (restricted) tree
    std::vector<std::string> matches;
    GetList( key, matches, section, subsection );

    // If there is no match, check for default
    if ( matches.size() == 0 ) {
      std::string defaultValue;
      flagStatus = CheckForDefault( defaultValue, key, section, subsection );

      // If there is a default, then test its value.
      // If it does not match, then re-set status
      if ( flagStatus == TRUE ) {
	if ( defaultValue != value ) {
	  flagStatus == FALSE;
	}
      }
    }

    // If there is a match, but it is not unique, call problem handler
    else if ( matches.size() > 1 ) {
      MultipleMatchHandler( key, section, subsection, matches.size() );
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


  // **************************************************************************
  //   Private Auxilliary Methods: Conversion Routines
  // **************************************************************************


  // =====================================
  //   Convert a DOMNode to a DOMElement
  // =====================================
  DOMElement* XMLParamHandler::Node2Elem( DOMNode *node ) {
    ENTER_FCN( "XMLParamHandler::Node2Elem" );
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
    ENTER_FCN( "XMLParamHandler::Node2Attr" );
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
  void XMLParamHandler::Node2String( DOMNode *textnode, std::string &value ){
    ENTER_FCN( "XMLParamHandler::Node2String" );
    if ( textnode->getNodeType() != DOMNode::TEXT_NODE ) {
      std::string errmsg;
      errmsg  = "Invalid conversion attempt from DOMNode to string!";
      errmsg += "\n Node is not of type TEXT_NODE\n Conversion of arbitrary";
      errmsg += " types to string is currently not supported!";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }
    else {
      // value = X2S(textnode->getNodeValue());
      value.assign(X2C(textnode->getNodeValue()));
    }
  }


  // **************************************************************************
  //   Private Auxilliary Methods: Search Routines
  // **************************************************************************


  // ========================================
  //   Recursive search method for elements
  // ========================================
  std::vector<DOMElement *>*
  XMLParamHandler::FindMatchingElements( std::vector<std::string> &keys,
					 DOMElement *treetop,
					 unsigned int curdepth ) {

    ENTER_FCN( "XMLParamHandler::FindMatchingElements" );

    // Perform consistency checks
    if ( keys.size() == 0 ) {
      Info->Error( "FindMatchingElements: Got zero keys!", __FILE__,
		   __LINE__);
    }
    if ( curdepth == 0 || curdepth > keys.size() ) {
      std::string errmsg = "FindMatchingElements: curdepth = ";
      errmsg += Info->GenStr( curdepth );
      errmsg += ", but should be in [1,";
      errmsg += Info->GenStr( keys.size() );
      errmsg += "]";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // If desired be verbose
    if ( beVerbose_ == true ) {
      std::cerr << " FindMatchingElements: curdepth = " << curdepth
		<< ", key[" << curdepth << "] = " << keys[curdepth-1]
		<< std::endl;
    }

    // Generate element vector
    std::vector<DOMElement *> *elemvec = new std::vector<DOMElement *>;

    // Generate list of all elements in subtree starting at treetop
    // that match the curdepth'th keyword
    DOMNodeList *list = treetop->getElementsByTagName(S2X(keys[curdepth-1]));

    // If we are on the lowest level allowed, convert the node list
    // to an element vector
    if ( curdepth == keys.size() ) {
      for ( unsigned int i = 0; i < list->getLength(); i++ ) {
	elemvec->push_back( Node2Elem( list->item(i) ) );
      }
    }

    // This is not the lowest level, so we descend and look for elements
    // in all subtrees spanned by the matches in our node list and append
    // the results to the element vector
    else {

      // auxilliary variable for intermediate results
      std::vector<DOMElement *> *tmpvec = NULL;

      for ( unsigned int i = 0; i < list->getLength(); i++ ) {

	// get results for each subtree
	tmpvec = FindMatchingElements( keys, Node2Elem(list->item(i)),
				       curdepth+1 );

	// append results to our element vector
	for ( unsigned int k = 0; k < tmpvec->size(); k++ ) {
	  elemvec->push_back( (*tmpvec)[k] );
	}

	// delete intermediate result to avoids memory leak
	delete tmpvec;
      }
      
    }

    // Return element vector
    return elemvec;
  }


  // ==========================================
  //   Recursive search method for attributes
  // ==========================================
  std::vector<DOMAttr*>*
  XMLParamHandler::FindMatchingAttributes( std::string attr_key,
					   std::vector<std::string> &keys,
					   DOMElement *treeTop,
					   std::vector<DOMElement*> *elemlist){

    ENTER_FCN( "XMLParamHandler::FindMatchingAttributes" );

    // Perform consistency check
    if ( keys.size() == 0 ) {
      std::string errmsg;
      errmsg  = "FindMatchingAttributes: Need at least one keyword! Got none!";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // Generate return vector
    std::vector<DOMAttr*> *attrvec = new std::vector<DOMAttr*>;

    // Clear element vector
    if ( elemlist != NULL ) {
      elemlist->clear();
    }

    // Find all matching elements
    std::vector<DOMElement *>* elements;
    elements = FindMatchingElements( keys, treeTop, 1 );

    // Be verbose, if demanded
    if ( beVerbose_ == true ) {
      std::cerr << " FindMatchingAttributes: Found " << elements->size()
		<< " matches for section '" << keys[0] << "'" << std::endl;
    }

    // Loop over all elements and examine their attributes
    DOMNamedNodeMap *attributes = NULL;
    DOMNode* aux_node = NULL;
    DOMAttr* match_attr = NULL;
    std::string dummy;
    for ( unsigned int i = 0; i < elements->size(); i++ ) {

      // Get attribute of elements
      attributes = (*elements)[i]->getAttributes();

      // Select attributes matching name and convert them
      aux_node = attributes->getNamedItem( S2X( attr_key ) );
      if ( aux_node != NULL ) {
	match_attr = Node2Attr( aux_node );
      }
      else {
	match_attr = NULL;
      }
      
      // If there is a match, append attribute to result vector
      // and also append element to element list
      if ( match_attr != NULL ) {
	attrvec->push_back( match_attr );
	if ( elemlist != NULL ) {
	  elemlist->push_back( (*elements)[i] );
	}
      }
    }

    // Finished
    return attrvec;

  }


  // ===============================================
  //   Return a list of strings matching a keyword
  // ===============================================
  void XMLParamHandler::FindAllMatches( const std::string key,
					std::vector<std::string> &list,
					const std::string section,
					const std::string subsection,
					DOMElement *treeTop,
					std::vector<DOMElement*> *elemlist ){

    ENTER_FCN( "XMLParamHandler::FindAllMatches" );

    // Report, what we are looking for
    if ( beVerbose_ == true ) {
      std::cerr << "\n FindAllMatches: Starting new search\n"
		<< "       key        = '" << key        << "'\n"
		<< "       section    = '" << section    << "'\n"
		<< "       subsection = '" << subsection << "'\n";
    }

    bool elem_match = false;
    bool attr_match = false;
    std::vector<DOMElement*>* elem_matches = NULL;
    std::vector<DOMAttr*   >* attr_matches = NULL;
    std::vector<DOMElement*> attr_elements;
    std::vector<std::string> keys;
    std::string value;

    // *************************************
    //   Part 1: Check of input parameters
    // *************************************

    // Make sure we got at least a key word
    if ( key == "" ) {
      Info->Error( "Error: Must specify at least one key word!", __FILE__,
		   __LINE__ );
    }

    // Check if vector is empty. If not issue a warning
    // and erase its entries, if this is desired
    if ( list.empty() != true ) {
      if ( beVerbose_ == true ) {
	std::string errmsg = "Warning input vector was not empty!\n";
	errmsg += "Contents have been erased!";
	Info->Warning( errmsg );
      }
      list.clear();
    }

    // Clear element vector
    if ( elemlist != NULL ) {
      elemlist->clear();
    }

    // **************************************************
    //   Part 2: Assume main key word is an element tag
    // **************************************************

    // Generate vector of keywords
    if ( section != "" ) {
      keys.push_back( section );
    }
    if ( subsection != "" ) {
      keys.push_back( subsection );
    }
    keys.push_back( key );

    // Find matching elements
    elem_matches = FindMatchingElements( keys, treeTop, 1 );

    // Check if there was a match
    if ( elem_matches->size() > 0 ) {
      elem_match = true;
    }

    // Be verbose, if demanded
    if ( beVerbose_ == true ) {
      std::cerr << " FindAllMatches: Found " << elem_matches->size()
		<< " matching elements for key '" << key << "'"<< std::endl;
    }

    // *****************************************************
    //   Part 3: Assume main key word is an attribute name
    // *****************************************************

    // Test whether a section was specified. If not, no attribute
    // search is possible
    if ( section != "" ) {

      // Generate vector of keywords
      keys.clear();
      keys.push_back( section );
      if ( subsection != "" ) {
	keys.push_back( subsection );
      }

      // Find matching attributes
      attr_matches = FindMatchingAttributes( key, keys, treeTop,
					     &attr_elements );

      // Be verbose, if demanded
      if ( beVerbose_ == true ) {
	std::cerr << " FindAllMatches: Found " << attr_matches->size()
		  << " matching attributes for key '" << key << "'"
		  << std::endl;
      }

      // Check if there was a match, if not call problem handler
      if ( attr_matches->size() > 0 ) {
	attr_match = true;
      }

    }

    // **********************************
    //   Part 4: Process search results
    // **********************************

    // See, if there were element and attribute matches
    if ( elem_match == true && attr_match == true ) {
      std::string errmsg = "Keyword '" + key + "' matches element(s) "
	+ "and attribute(s)!";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // Convert element values to strings and append to list
    // If desired also save found elements
    else if ( elem_match == true ) {
      for ( unsigned int i = 0; i < elem_matches->size(); i++ ) {
	value.assign( GetElementValue( (*elem_matches)[i] ) );
	list.push_back( value );
	if ( elemlist != NULL ) {
	  elemlist->push_back( (*elem_matches)[i] );
	}
      }
    }

    // Convert element values to strings and append to list
    // If desired also save found elements
    else if ( attr_match == true ) {
      for ( unsigned int i = 0; i < attr_matches->size(); i++ ) {
	value.assign( X2C( (*attr_matches)[i]->getValue() ) );
	list.push_back( value );
	if ( elemlist != NULL ) {
	  elemlist->push_back( attr_elements[i] );
	}
      }
    }


    // ******************
    //   Part 5: Report
    // ******************
    if ( LOGTOINFO ) {
      std::string msg = "Section '" + section + "' -> Subsection '" +
	subsection + "' -> Keyword '" + key + "'";
      Info->PrintF( "", "%s", msg.c_str() );
      for ( unsigned int k = 0; k < list.size(); k++ ) {
	msg = "Value = '" + list[k] + "'";
	Info->PrintF( "", "%s", msg.c_str() );
      }
    }

    // *******************
    //   Part 6: Cleanup
    // *******************
    delete elem_matches;
    delete attr_matches;

  }


  // **************************************************************************
  //   Private Auxilliary Methods: Treatment of errors and defaults
  // **************************************************************************


  // ==================================
  //   Treat case of multiple matches
  // ==================================
  void XMLParamHandler::MultipleMatchHandler( const std::string key,
					      const std::string section,
					      const std::string subsection,
					      const unsigned int nmatches ) {

    ENTER_FCN( "XMLParamHandler::MultipleMatchHandler" );

    std::string errmsg;

    // Make sure that we have a non-unique match
    if ( nmatches > 1 ) {
      errmsg += "Error: Match for keyword '" + key + "' is not unique";
      if ( section != "" && subsection == "" ) {
	errmsg += " within sections '" + section + "'";
      }
      if ( subsection != "" ) {
	errmsg += " within subsections '" + subsection;
	errmsg += " of sections '" + section + "'";
      }
      errmsg += '\n';
    }

    // This should not happen!
    else {
      errmsg  = "Internal Error: MultipleMatchHandler called with nmatches = ";
      errmsg += Info->GenStr(nmatches);
      errmsg += "!\n This should not have happened!\n";
    }

    // And now terminate
    Info->Error( errmsg, __FILE__, __LINE__ );

  }


  // =====================================
  //   Treat case of no match for string
  // =====================================
  void XMLParamHandler::NoMatchHandler( std::string &value,
					const std::string key,
					const std::string section,
					const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::NoMatchHandler" );

    // Test, whether a default value is specified for the parameter
    std::string defaultValue;
    Boolean defaultExists = CheckForDefault( defaultValue, key, section,
					     subsection );

    // If default exist, return it
    if( defaultExists == TRUE ) {
      value = defaultValue;
    }

    // No match and no default value, so cry out!
    else {
      NoMatchErrorReporter( key, section, subsection );
    }

  }


  // ======================================
  //   Treat case of no match for Integer
  // ======================================
  void XMLParamHandler::NoMatchHandler( Integer &value,
					const std::string key,
					const std::string section,
					const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::NoMatchHandler" );

    // Test, whether a default value is specified for the parameter
    std::string defaultValue;
    Boolean defaultExists = CheckForDefault( defaultValue, key, section,
					     subsection );

    // If default exist, convert it to Integer
    if( defaultExists == TRUE ) {
      value = String2Integer( defaultValue );
    }

    // No match and no default value, so cry out!
    else {
      NoMatchErrorReporter( key, section, subsection );
    }

  }


  // =====================================
  //   Treat case of no match for Double
  // =====================================
  void XMLParamHandler::NoMatchHandler( Double &value, const std::string key,
					const std::string section,
					const std::string subsection ) {

    ENTER_FCN( "XMLParamHandler::NoMatchHandler" );

    // Test, whether a default value is specified for the parameter
    std::string defaultValue;
    Boolean defaultExists = CheckForDefault( defaultValue, key, section,
					     subsection );

    // If default exist, convert it to Integer
    if( defaultExists == TRUE ) {
      value = String2Double( defaultValue );
    }

    // No match and no default value, so cry out!
    else {
      NoMatchErrorReporter( key, section, subsection );
    }

  }


  // =================================
  //   Report that there is no match
  // =================================
  void XMLParamHandler::NoMatchErrorReporter( const std::string key,
					      const std::string section,
					      const std::string subsection ){

    ENTER_FCN( "XMLParamHandler::NoMatchErrorReporter" );

    std::string errmsg;

    errmsg += "Error: No match and no default found for keyword '" + key
      + "'";
    if ( section != "" && subsection == "" ) {
      errmsg += " within sections '" + section + "'";
    }
    if ( subsection != "" ) {
      errmsg += " within subsections '" + subsection + "'";
      errmsg += " of sections '" + section + "'";
    }
    errmsg += '\n';
    Info->Error( errmsg, __FILE__, __LINE__ );
  }


  // =================================
  //   Check for default parameter
  // =================================
  Boolean XMLParamHandler::CheckForDefault( std::string &defaultValue,
					    const std::string key,
					    const std::string section,
					    const std::string subsection ) {

    Boolean defaultFound = FALSE;

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

    // Find all elements/values matching keyword in (restricted) defaults
    // tree
    std::vector<std::string> matches;
    FindAllMatches( key, matches, section, subsection, rootElemDefaults_ );

    // If there was no unique match, call problem handler
    if ( matches.size() > 1 ) {
      MultipleMatchHandler( key, section, subsection, matches.size() );
    }
  
    // Check, if a default was found
    if ( matches.size() == 1 ) {
      defaultFound = TRUE;
      defaultValue = matches[0];
    }

    // If no default was found, test, whether this is one of the special
    // cases not implemented so far
    else {

      if ( key == "mesh_library" ) {
	defaultFound = TRUE;
	defaultValue = "cfsgrid";
      }

      if( key == "preStressVal" ) {
	defaultFound = TRUE;
	defaultValue = "0";
      }

      if( key == "effMass" ) {
	defaultFound = TRUE;
	defaultValue = "no";
      }

    }

    // Tell what we found
    if ( beVerbose_ == true ) {
      std::string msg = "CheckForDefault: Default for parameter '" + key;
      msg += "' is '" + defaultValue + "'";
      Info->Warning( msg );
    }

    return defaultFound;

  };


  // **************************************************************************
  //   Private Auxilliary Methods: Miscellaneous
  // **************************************************************************


  // ===========================
  //   Get value of an element
  // ===========================
  char* XMLParamHandler::GetElementValue( DOMElement *elem ) {

    ENTER_FCN( "XMLParamHandler::GetElementValue" );

    // Obtain child nodes of this element
    DOMNodeList *children = elem->getChildNodes();

    // Complain, if there is not exactly one child
    if ( children->getLength() != 1 ) {
      std::string errmsg;
      if ( children->getLength() == 0 ) {
	errmsg = "GetElementValue: Encountered element without child!";
      }
      else {
	errmsg ="GetElementValue: Encountered element with multiple ";
      }
      errmsg += "children!\nGetElementValue: Element tag is '";
      errmsg += X2C(elem->getNodeName());
      errmsg += "'";
      Info->Error( errmsg, __FILE__, __LINE__ );
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

    ENTER_FCN( "XMLParamHandler::GetElementAttribute" );

    bool hasAttr;
    DOMNamedNodeMap *attributes = NULL;
    DOMNode* aux_node = NULL;
    DOMAttr* aux_attr = NULL;

    // Get attribute of elements
    attributes = element->getAttributes();

    // Select attributes matching name and convert them
    aux_node = attributes->getNamedItem( S2X( keyword ) );

    // Check, if element has specified attribute
    if ( aux_node != NULL ) {
      hasAttr = true;
    }
    else {
      hasAttr = false;
    }

    // If element has specified attribute obtain its value
    if ( hasAttr ) {
      aux_attr = Node2Attr(aux_node);
      attrVal.assign( X2C( aux_attr->getValue() ) );
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

    // String for assembling error messages
    std::string errmsg;

    // Create the parser, force validation using full schema checking with
    // namespaces and tell it to drop irrelevant white space
    (*parser) = new XercesDOMParser();
    (*parser)->setValidationScheme(XercesDOMParser::Val_Always);
    (*parser)->setDoNamespaces(true);
    (*parser)->setDoSchema(true);
    (*parser)->setValidationSchemaFullChecking(true);
    (*parser)->setIncludeIgnorableWhitespace(false);

    // We may separate the schema file from the instance file
    std::string cfsSchema = "http://www.cfs++.org ";
    cfsSchema += CVSEXTERNAL;
    cfsSchema += "/CFS++XML/CFS.xsd";
    (*parser)->setExternalSchemaLocation( cfsSchema.c_str() );

    // Have not yet understood what an entity reference node is, but it seems
    // we do not need them
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
      errmsg  = "The following XML error was encountered during parsing:\n";
      errmsg += X2C( event.getMessage() );
      Info->Error( errmsg, __FILE__, __LINE__ );
    }
    catch( const DOMException &event ) {
      const unsigned int maxChars = 2047;
      XMLCh errText[maxChars + 1];

      errmsg = "DOM Error during parsing! DOMException code is:\n"
	+ Info->GenStr( event.code ) + '\n';

      if( DOMImplementation::loadDOMExceptionMsg(event.code, errText,
						 maxChars) ) {
	errmsg += "Message is: ";
	errmsg += X2C(errText);
      }

      Info->Error( errmsg, __FILE__, __LINE__ );
    }
    catch(...) {
      errmsg  = "An unknown error occurred during parsing!\n";
      errmsg += "All I can say is that it was neither an XMLException"; 
      errmsg += " nor a DOMException.";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // Obtain and validate root element of document tree
    DOMDocument *doc = (*parser)->getDocument();
    DOMNodeList *children = doc->getChildNodes();
    if ( children->getLength() != 1 ) {
      errmsg  = "Document root has more than one child!\n";
      errmsg += "We are in real trouble here!";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }
    if ( children->item(0)->getNodeType() != DOMNode::ELEMENT_NODE ) {
      errmsg  = "Root node is not of type DOMNode::ELEMENT_NODE!\n";
      errmsg += "We are in real trouble here!";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }
    rootElem = (DOMElement *)(children->item(0)); 

    // We are finished
    return rootElem;
  }

}

#endif
