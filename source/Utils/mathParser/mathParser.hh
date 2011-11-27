#ifndef MATHPARSER_HH
#define MATHPARSER_HH

#include <map>
#include <list>
#include <set>
#include <string>
#include <boost/signals.hpp>
#include "muParser.h"
#include "Utils/StdVector.hh"
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "General/Environment.hh"

namespace CoupledField {

  template <class TYPE> class Vector;

  //! Forward class declarations
  class CoordSystem;
  
  //! Handles mathematical parser for different contexts
  class MathParser{

  public:

    //! Datatype for identifying a specific mathParser
    typedef UInt HandleType;

    //! Abbreviation for global handle
    enum {GLOB_HANDLER=0};
    
    //! Constructor
    MathParser();
    
    //! Destructor
    ~MathParser();

    //! Get new parser handle

    //! Returns a handle for a new parser. This handle is needed for all
    //! subsequent communication with this parser.
    //! \param setDefaults If true, the parser will set default values
    //!                    for the most common variables (f,t,x,y,z), which
    //!                    is useful if for an expression the global values
    //!                    are not yet set, but will be defined at a later point.
    //!                    If a default variable is already set in the GLOBAL
    //!                    parser instance, this value will be taken. 
    //! \return New MathParser handle
    
    HandleType GetNewHandle( bool setDefaults = false);

    //! Free the memory and parameters associated with a certain handle

    //! This method frees the internal parser, its memory and variables 
    //! associated with this handle.
    //! \param handle Parser handle identifying an internal parser(context)
    //!               to be freed.
    void ReleaseHandle( HandleType handle );

    //! Pass the expression to be evaluated to the parser
    
    //! This methods passes an expression to be evaluated to the parser
    //! itself. This will trigger the syntactical analysis of the 
    //! expression. In order to evaluate the expression, a successive call
    //! to Eval() has to be performed.
    //! The expression can be composed of several subexrepssions, 
    //! delimited by ','. In this case, all the subexpressions can be 
    //! retrieved using either EvalVector() or EvalMatrix().
    //! \param handle MathParser handle for identifying specific parser
    //! \param expr Expression to be parsed
    void SetExpr( HandleType handle, const std::string &expr );
    
    //! Retrieve the expression of the related parser

    //! This methods returns the parser expression of the related parser
    //! instance.
    //! \param handle MathParser handle for identifying specific parser
    //! \return Expression of the parsewr
    std::string GetExpr( HandleType handle );

    //! Query if current expression is constant
    
    //! Returns if the current expression of the related parser instance
    //! is constant (not depending on any variables)
    //! \param handle MathParser handle for identifying specific parser
    //! \return Flag is expression is constant
    bool IsExprConstant( HandleType handle );
    
    //! Query if current expression uses specified variable \a var
    
    //! Returns, is the current expression of the related parser instance
    //! depends on the variable .
    //! \param handle MathParser handle for identifying specific parser
    //! \param var Variable name the expression is questioned for
    //! \return Flag if expression contains variable
    bool IsExprVariable( HandleType handle, const std::string& var );
    
    //! Get all variables the expression of the parser instance depends on
    
    //! This method returns all the variables, which are used in the current
    //! expression of the specified parser instance.
    //! \param handle MathParser handle for identifying specific parser
    //! \param varNames List of variables the expression depends on
    void GetExprVars( HandleType handle, StdVector<std::string>& varNames );
    
    
    //! Return number of expressions set
    
    //! This method returns the number of expression set for the parser
    //! denoted by #handle. 
    UInt GetNumExprs( HandleType handle );
    
    //! Evaluate mathematical expression previously set by SetExpr()

    //! This method evaluates the expression previously set by SetExpr().
    //! If several, comma separated expressions are set, only the last
    //! one is taken.
    //! \param handle MathParser handle for identifying specific parser
    Double Eval( HandleType handle );
    
    //! Evaluate list of mathematical expression, set by SetExpr()
    
    //! This method evaluated the expression previously set by SetExpr()
    //! and stores the values in the vector.
    //! \param handle MathParser handle for identifying specific parser
    //! \param vec Vector to be filled
    void EvalVector( HandleType handle, Vector<Double>& vec );
    
    //! Evaluate expressions to matrix, previously set by SetExpr()
    
    //! This method evaluates the expression previously set by SetExpr()
    //! and stores the values in the matrix.
    //! The function tries to be smart about the size of the matrix
    //! as follows:
    //! - If #numRows and / or #numCols are given, the matrix gets resized
    //!   to the number of rows / cols given. If the number of entries does
    //!   not match the final size of the matrix, an exception is thrown.
    //! - If neither #numRows nor #numCols is prescribed, the initial size
    //!   of the matrix is taken. If the number of entries does not match 
    //!   the matrix size, an exception is thrown.
    //! If neither #numRows nor #numCols are given and the matrix has zero
    //! size, an exception is thrown as well.
    //! \param handle MathParser handle for identifying specific parser
    //! \param matrix A matrix where the values get stored into
    //! \param numRows Number of rows, if 
    void EvalMatrix( HandleType handle, Matrix<Double>& matrix,
                     UInt numRows = 0, UInt numCols = 0 ); 
    
    //! Dump all parser instances, their variables and expression
    void Dump( std::ostream& os );

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
    //! \param handle MathParser handle for identifying specific parser
    //! \param varName Name of variable to be set
    //! \param val Value of variable
    void SetValue( HandleType handle,
                   const std::string &varName,
                   Double val );
    
    //! Special Set function for coordinate values

    //! This function sets the values of the coordinate components
    //! of the given vector within the given parser.
    //! For a cylindrical coordinate system for example, afterwards the
    //! following expression is possible:
    //! \verbatim
    //!  z * 10 + sin( phi )
    //! \endverbatim
    //! \param handle MathParser handle for identifying specific parser
    //! \param coosy Local coordinate system, for which the coordinate components
    //!              are to be registered within the given parser
    //! \param globCoord Global coordinates (x,y,z) of the given point
    void SetCoordinates( HandleType handle,
                         const CoordSystem &coosy,
                         const Vector<Double> &globCoord );
    //@}
    
    // =======================================================================
    //  CALLBACK FUNCTIONALITY
    // =======================================================================
    //@{
        //! \name Callback functionality of mathParser
        
    //! Define signal type of MathParser
    typedef boost::signal<void()> MathParserSignal;
    
    //! Register callback function for change of value of expression
    boost::signals::connection 
    AddExpChangeCallBack( const MathParserSignal::slot_function_type
                          &subscriber,
                          HandleType handle );
    //@}
    
  protected:

    //! Typedef for variable pool (for one handle )
    typedef std::map<std::string, Double>  VarPool;
    
    //! Typedef for mapping handle to variable pool
    typedef std::map<HandleType, VarPool> PoolMap;
    
    //! Typedef for parser pool
    typedef std::map<HandleType, mu::Parser> ParserMap;
    
    //! Initialize parser
    void InitParser( mu::Parser &parser, VarPool &pool, 
                     bool isGlobal,
                     bool setDefaults );

    //! Get parser by handle
    mu::Parser & GetParser( HandleType handle );
    
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
    
    //! Map for every global variable parser instance which uses this variable
    std::map<std::string, std::set<HandleType> > globVarsInUse_;
    
    //! Set with used variables in each expression 
    std::map<HandleType, std::set<std::string > > varsInUse_;

    //! Set with currently active handles
    std::set<HandleType> activeHandles_;

    //! Memory for implicitly allocated variables
    static std::list<Double> dynamicPool_;
    
    // =======================================================================
    //  SIGNAL HANDLING STUFF
    // =======================================================================
    
    //! Type definition for shared pointer to MathParserSignal
    typedef shared_ptr<MathParserSignal> PtSig;
    
    //! Callback signals if value of expression changes
    std::map<HandleType, PtSig > exprChangeSignal_;
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
