// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef MATHPARSER_HH
#define MATHPARSER_HH

#include <map>
#include <list>
#include "muParser.h"
#include "Utils/StdVector.hh"
#include "General/environment.hh"


namespace CoupledField {

  //! Forward class declarations
  class CoordSystem;
  
  //! Handless mathematical parser for different contexts
  class MathParser{

  public:

    //! Datatype for identifying a specific mathParser
    typedef UInt HandleType;

    //! Abbreviation for global handler
    enum {GLOB_HANDLER=0};
    
    //! Constructor
    MathParser();
    
    //! Destructor
    ~MathParser();

    //! Get new parser handle

    //! Returns a handle for a new parser. This handle is needed for all
    //! subsequent communication with this parser.
    //! \return New MathParser handler
    HandleType GetNewHandle();

    //! Free the memory and parameters associated with a certain handle

    //! This method frees the internal parser, its memory and variables 
    //! associated with this handle.
    //! \param handle Parser handle identifying an internal parser(context)
    //!               bo be freed.
    void ReleaseHandle( HandleType handle );

    //! Pass the expression to be evaluated to the parser
    
    //! This methods passes an expression to be evaluated to the parser
    //! itself. This will trigger the syntactical analysis of the 
    //! expression. In order to evaluate the expression, a successive call
    //! to Eval() has to be performed
    void SetExpr( HandleType handler, const std::string &expr );

    //! Evaluate mathematical expression previously set by SetExpr()

    //! This method evaluates the expression previously set by SetExpr()
    //! \param handler MathParser handler for identifying specific parser
    Double Eval( HandleType handler );

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
    void SetValue( HandleType handler,
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
    void SetCoordinates( HandleType handler,
                         const CoordSystem &coosy,
                         const Vector<Double> &globCoord );
    //@}

  protected:

    //! Typedef for variable pool (for one handler )
    typedef std::map<std::string, Double>  VarPool;
    
    //! Typedef for mapping handler to variable pool
    typedef std::map<HandleType, VarPool> PoolMap;
    
    //! Typedef for parser pool
    typedef std::map<HandleType, mu::Parser> ParserMap;
    
    //! Initialize parser
    void InitParser( mu::Parser &parser, bool isGlobal );

    //! Get parser by handler
    mu::Parser & GetParser( HandleType handler );
    
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

    // =======================================================================
    //  CFS SPECIFIC ADDITIONAL FUNCTIONS
    // =======================================================================
    
    //! Get local coordinate component of coordinate system with given id (3d)
    static Double LocCoord3D( const char * coordSysId, 
                              Double dof,
                              Double x, Double y, Double z );
    
    //! Get local coordinate component of coordinate system with given id (2d)
    static Double LocCoord2D( const char * coordSysId, 
                              Double dof,
                              Double x, Double y);

        
    //! Pool of math parsers
    ParserMap parsers_;
    
    //! Pool of variables
    PoolMap pools_;

    //! Set with currently active handles
    std::set<HandleType> activeHandles_;

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
  //! The class implements the following concepts:
  //! - There can exist several parser objects within a program with different
  //!   known variables and functions, which are all
  //!   administered by the MathParser-class. To identify a specific parser, a
  //!   MathParserHandle is used, which uniquely identifies a certain parser
  //!   context
  //! - There exist two sets of administered parser (identified by their 
  //!   Handle):
  //!   One consists only of the global parser
  //!   ( defined by the Handle GLOB_HANDLER ) and the other only of 
  //!   normal/local ones, which can
  //!   be created on demand. The difference between both sets is, that all
  //!   variables set in the global parser will automatically also be promoted
  //!   to all the other parser objects but not vice-versa. A global variable
  //!   is for example the current (global) timestep t, whereas to coordinates
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
