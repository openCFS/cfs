// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_xerces.hh>

#ifdef USE_XERCES

#ifndef XERCES_HH_
#define XERCES_HH_

#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>

#include "DataInOut/ParamHandling/ParamNode.hh"


namespace CoupledField {

  class ParamNode;

  /** This class fills a ParamNode element with an XML file content using the DOM parser from xerces. 
   * <p>This might be a good place for any XML stuff besides only filling a ParamNode</p>
   * <p>It might be a good idea, to release this object soon to free the system ressources. Note, that
   * the Xerces object is not needed to keep the ParamNode createted by CreateParamNodeInstance().</p>
   * @see PtrParamNode*/
  class Xerces 
  {

    public:
      /** This sets the internal variables but dos not parse the file yet.
       * @param schema if given a validation is done. */
      Xerces(const std::string& file, const std::string& schema = "");

      /** Shuts down all xerces stuff */ 
      ~Xerces();
      
      /** creates a param node instance filled with the xml content from the constuctor.
       * Easy to modify for using an arbitrary DOMNode* data source
       * @return the caller has to delete the tree by itself, there is no reference within the Xerces class.  */ 
      PtrParamNode CreateParamNodeInstance();

    private:
      /** This actually parses the file and sets root_ and the other variables */
      void Parse();

      /** recursive implementation of Fill */
      void Fill(xercesc::DOMNode* node, PtrParamNode out);

      /** The result of the constructor */    
      xercesc::DOMNode* root_;
      
      /** actually not needed any more after the constructor but deletet only in the destructor. */      
      xercesc::XercesDOMParser* parser_;

      /** The filename we are based on */
      std::string file_;
      
      /** The schema filename is empty if not given */
      std::string schema_; 

      /** This is an implementation of the xerces error handler, the methods are called if parsing
       * stumbles over the file/schema. The implementation is simply a logging */
      class EventHandler: public xercesc::ErrorHandler
      {
        public:
          /** @param xerces the object we handle the events for -> required for nice outputs */
          EventHandler(const Xerces* xerces);
          
          /** logs */
          void warning( const xercesc::SAXParseException &event );
      
          /** throws an Exception */
          void error( const xercesc::SAXParseException &event );
      
          /** throws an Exception */
          void fatalError( const xercesc::SAXParseException &event );
          
          /** not implemented but needs to be there as it is abstract in the base class */
          void resetErrors();
          
        private:
          /** This is the object we handle the events */
          const Xerces* xerces_;  
          
      }; // end of EventHandler     
      
  }; // end of Xerces


} // end of namespace


#endif // XERCES_HH_

#endif // USE_XERCES
