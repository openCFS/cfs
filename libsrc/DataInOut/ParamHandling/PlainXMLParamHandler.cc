
//#include "General/environment.hh"
//#include "DataInOut/WriteInfo.hh"

#include "BaseParamHandler.hh"
#include "PlainXMLParamHandler.hh"


namespace CoupledField {


  // **************************************************************************
  //   Constructors and Destructors
  // **************************************************************************

  // ===============
  //   Constructor
  // ===============
  PlainXMLParamHandler::PlainXMLParamHandler( char *fname ) {}


  // ==============
  //   Destructor
  // ==============
  // PlainXMLParamHandler::~PlainXMLParamHandler() {}


  // **************************************************************************
  //   Public Methods
  // **************************************************************************


  // ================================================================
  //   Return as string the value of a parameter matching a keyword
  // ================================================================
  void PlainXMLParamHandler::Get( const std::string key, std::string &value,
				  const std::string section,
				  const std::string subsection ) {}


  // =================================================================
  //   Return as integer the value of a parameter matching a keyword
  // =================================================================
  void PlainXMLParamHandler::Get( const std::string key, int &value,
				  const std::string section,
				  const std::string subsection ) {}


  // ================================================================
  //   Return as double the value of a parameter matching a keyword
  // ================================================================
  void PlainXMLParamHandler::Get( const std::string key, double &value,
				  const std::string section,
				  const std::string subsection ) {}


  // ====================================================
  //   Return list of strings values matching a keyword
  // ====================================================
  void PlainXMLParamHandler::GetList( const std::string key,
				      std::vector<std::string> &list,
				      const std::string section,
				      const std::string subsection ) {}


  // ===================================================
  //   Return list of Double values matching a keyword
  // ===================================================
  void PlainXMLParamHandler::GetList( const std::string key,
				      std::vector<Double> &list,
				      const std::string section,
				      const std::string subsection ) {}

  // =====================================
  //   Return a list of the defined PDEs
  // =====================================
  void PlainXMLParamHandler::GetPDEList( std::vector<std::string> &list ) {}


  // ========================================
  //   Query on/off status of a flag/switch
  // ========================================
  Boolean PlainXMLParamHandler::IsSet( const std::string key,
				       const std::string section,
				       const std::string subsection )
  {
    return FALSE;
  }


  // =================================================
  //   Query whether a parameter has a certain value
  // =================================================
  Boolean PlainXMLParamHandler::HasValue( const std::string key,
					  const std::string value,
					  const std::string section,
					  const std::string subsection )
  {
    return FALSE;
  }

}
