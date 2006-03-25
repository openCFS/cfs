#ifndef MATHPARSER_HH
#define MATHPARSER_HH

#include <map>
#include <list>
#include "muParser.h"
#include "StdVector.hh"
#include "General/environment.hh"


namespace CoupledField {

  //! Forward class declarations
  class CoordSystem;
  
  //! Handless mathematical parser for different contexts
  class MathParser{

  public:

    //! Datatype for identifying a specific mathParser
    typedef UInt HandlerType;

    //! Abbreviation for global handler
    enum {GLOB_HANDLER=0};
    
    //! Constructor
    MathParser();
    
    //! Destructor
    ~MathParser();

    //! Get new parser handler

    //! Returns a handler for a new parser. This handler is needed for all
    //! subsequent communication with this parser.
    //! \return New MathParser handler
    HandlerType GetNewHandler();

    //! Evaluate mathematical expression

    //! This method evaluates the given expression using in the specified
    //! parser.
    //! \param handler MathParser handler for identifying specific parser
    //! \param expr Mathematical expression to be evaluated
    Double Eval( HandlerType handler, const std::string &expr );

    // =======================================================================
    //  SET METHODS
    // =======================================================================

    //@{
    //! \name Set Methods For Additional Variables

    //! Introduce a variable within the parser with given value

    //! This methods registers the given variable name within the specified 
    //! parser and assigns it the given value. If the function is called for
    //! first time, memory is allocated,the variable name is registered in
    //! the parser and the value ist set. All subsequent calls just update the
    //! variable value.
    //! \param handler MathParser handler for identifying specific parser
    //! \param varName Name of variable to be set
    //! \param val Value of variable
    void SetValue( HandlerType handler,
                   const std::string &varName,
                   Double val );
    
    //! Special Set function for coordinate values

    //! This function sets the values of the coordinate components
    //! of the given vector within the given parser.
    //! For a cylindrical coordinate system for example, afterwards the
    //! following expression is possible:
    //! \verbatim
    //!  uax * 10 + sin( uphi )
    //! \endverbatim
    //! \param handler MathParser handler for identifying specific parser
    //! \param coosy Local coordinate system, for which the coordinate components
    //!              are to be registered within the given parser
    //! \param globCoord Global coordinates (x,y,z) of the given point
    void SetCoordinates( HandlerType handler,
                         const CoordSystem &coosy,
                         const Vector<Double> &globCoord );
    //@}

  protected:

    //! Typedef for variable pool (for one handler )
    typedef std::map<std::string, Double>  VarPool;
    
    //! Typedef for mapping handler to variable pool
    typedef std::map<HandlerType, VarPool> PoolMap;
    
    //! Typedef for parser pool
    typedef std::map<HandlerType, mu::Parser> ParserMap;
    
    //! Initialize parser
    void InitParser( mu::Parser &parser, Boolean isGlobal );

    //! Get parser by handler
    mu::Parser & GetParser( HandlerType handler );
    
    //! Factory function for adding variables
    static Double * AddVariable( const char *varName );

    // =======================================================================
    //  HELPER METHODS
    // =======================================================================   
    //@{
    //! \name Static Aliases For Operators <,<=, >, >=
    static Double Op_gt( Double a, Double b ) { return a > b; }
    static Double Op_ge( Double a, Double b ) { return a >= b; }
    static Double Op_le( Double a, Double b ) { return a <= b; }
    static Double Op_lt( Double a, Double b ) { return a < b; }
    //@}
        
    //! Pool of math parsers
    ParserMap parsers_;
    
    //! Pool of variables
    PoolMap pools_;

    //! Memory for implicitly allocated variables
    static std::list<Double> dynamicPool_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MathParser
  //! 
  //! \purpose 
  //! This class allows one to evaluate mathematical expressions, given in 
  //! string representation, in different context, i.e. in different 
  //! (internally administered) contexts, with different known variables. 
  //! Therefore is utilizes the open-source library muParser,
  //! a highly efficient parser library, which internally stores
  //! expressions in bytecode format. Additionally one can define own 
  //! variables and functions within the parser object.
  //!
  //! The class implements the following concepts
  //! - There can exist several parser objects within a program with different
  //!   known variables and functions, which are all
  //!   administered by the MathParser-class. To identify a specific parser, a
  //!   MathParserHandle is used.
  //! - There exist two set of administered parser (identified by their 
  //!   Handler):
  //!   One consists only of the global parser
  //!   ( defined by the Handler GLOB_HANDLER ) and the other only of 
  //!   normal/local ones, which can
  //!   be created on demand. The difference between both sets is, that all
  //!   variables set in the global parser will automatically also be promoted
  //!   to all the other parser object,s but not vice-versa. A global variable
  //!   is for example the current (global) timestep, whereas to coordinates
  //!   of the current boundary node are defined only in a local parser 
  //!   context within a specific PDE.
  //! 
  //! \collab 
  //! This class relies on the muParser library. Also it provides convenient
  //! functions for automatically setting the local representation of a vector
  //! associated to a local coordinate system.
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

}

#endif
