#ifndef MATHPARSER_HH
#define MATHPARSER_HH

#include <map>
#include <list>
#include <set>
#include <string>
#include <boost/signals2.hpp>
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

    //! The handle type for the muParser is simply unsigned int

    // This is the same as MathParser_GLOB_HANDLER in EnvironmentTypes.hh
    enum { GLOB_HANDLER = 0 };

    //! Constructor
    MathParser();
    
    //! Destructor
    virtual ~MathParser();

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
    
    virtual unsigned int GetNewHandle( bool setDefaults = false);

    //! Free the memory and parameters associated with a certain handle

    //! This method frees the internal parser, its memory and variables 
    //! associated with this handle.
    //! \param handle Parser handle identifying an internal parser(context)
    //!               to be freed.
    virtual void ReleaseHandle( unsigned int handle );

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
    virtual void SetExpr( unsigned int handle, const std::string &expr );
    
    //! Retrieve the expression of the related parser

    //! This methods returns the parser expression of the related parser
    //! instance.
    //! \param handle MathParser handle for identifying specific parser
    //! \return Expression of the parser
    std::string GetExpr( unsigned int handle );

    //! Query if current expression is constant
    
    //! Returns if the current expression of the related parser instance
    //! is constant (not depending on any variables)
    //! \param handle MathParser handle for identifying specific parser
    //! \return Flag is expression is constant
    bool IsExprConstant( unsigned int handle );
    
    //! Query if current expression uses specified variable \a var
    
    //! Returns, is the current expression of the related parser instance
    //! depends on the variable .
    //! \param handle MathParser handle for identifying specific parser
    //! \param var Variable name the expression is questioned for
    //! \return Flag if expression contains variable
    bool IsExprVariable( unsigned int handle, const std::string& var );
    
    //! Get all variables the expression of the parser instance depends on
    
    //! This method returns all the variables, which are used in the current
    //! expression of the specified parser instance.
    //! \param handle MathParser handle for identifying specific parser
    //! \param varNames List of variables the expression depends on
    void GetExprVars( unsigned int handle, StdVector<std::string>& varNames );

    //! Get variable of expression of the parser instance depends on

    //! This method returns value of a certain variable.
    //! \param handle MathParser handle for identifying specific parser
    //! \param varName variable to het value of
    Double GetExprVars( unsigned int handle, std::string varName );

    //! Get pointer to variable value for fast repeated access

    //! This method returns a pointer to the stored value of a variable,
    //! allowing fast repeated access without map lookups.
    //! \param handle MathParser handle for identifying specific parser
    //! \param varName Name of variable to get pointer to
    //! \return Pointer to the stored Double value, or nullptr if not found
    //! \note The pointer remains valid as long as the variable exists
    Double* GetValuePtr( unsigned int handle, const std::string& varName );

    //! Return number of expressions set
    
    //! This method returns the number of expression set for the parser
    //! denoted by #handle. 
    UInt GetNumExprs( unsigned int handle );
    
    //! Evaluate mathematical expression previously set by SetExpr()

    //! This method evaluates the expression previously set by SetExpr().
    //! If several, comma separated expressions are set, only the last
    //! one is taken.
    //! \param handle MathParser handle for identifying specific parser
    virtual Double Eval( unsigned int handle );
    
    //! Evaluate list of mathematical expression, set by SetExpr()
    
    //! This method evaluated the expression previously set by SetExpr()
    //! and stores the values in the vector.
    //! \param handle MathParser handle for identifying specific parser
    //! \param vec Vector to be filled
    virtual void EvalVector( unsigned int handle, Vector<Double>& vec );
    
    //! Evaluate list of mathematical expression, set by SetExpr()

    //! This method evaluated the expression previously set by SetExpr()
    //! and computes the divergence of it
    //! \param handle MathParser handle for identifying specific parser
    //! \param divergence result of computation
    virtual void EvalDivVector( unsigned int handle, Double& divergence );

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
    virtual void EvalMatrix( unsigned int handle, Matrix<Double>& matrix,
                     UInt numRows = 0, UInt numCols = 0 ); 
    
    /** gives all variables and their values for a handle. This is part of the information given by Dump().
     * @param handle most interesting probably for MathParser_GLOB_HANDLER. Nothing done for invalid handle
     * It would be nice to have a noexept but an icc bug forbids it :( https://software.intel.com/en-us/node/629136 */
    StdVector<std::pair<std::string, double> > GetRegisteredValues(unsigned int handle) const; // icc bug prevents noexcept

    /** for debug purpose: comma separated list of GetRegisteredValues() content */
    std::string ToString(unsigned int handle) const;

    /** creates a comma separated list of registered variables */
    std::string GetRegisteredVariables(unsigned int handle) const; // icc bug prevents noexcept

    /** reports the variables and value
     * @param handle nothing done for invalid handle */
    void ToInfo(PtrParamNode pn, unsigned int handle) const;

    //! Dump all parser instances, their variables and expression
    virtual void Dump( std::ostream& os );

    // =======================================================================
    //  SET METHODS
    // =======================================================================

    //@{
    //! \name Set Methods For Additional Variables

    //! Introduce a variable within the parser with given value

    //! This methods registers the given variable name within the specified 
    //! parser and assigns it the given value. If the function is called for
    //! first time, memory is allocated,the variable name is registered in
    //! the parser and the value is set. All subsequent calls just update the
    //! variable value.
    //! \param handle MathParser handle for identifying specific parser
    //! \param varName Name of variable to be set
    //! \param val Value of variable
    virtual void SetValue( unsigned int handle,
                   const std::string &varName,
                   Double val );
    
    //! Register external variable with given parser
    
    //! This method allows to register external variables, which are not handled
    //! by the internal variable pool of the MathParser. This also means, that
    //! no signal is fired, if the value changes, as it is the case with the method
    //! MathParser::SetValue(). On the other hand, this method has less overhead
    //! and should perform faster.
    //! \param handle MathParser handle for identifying specific parser
    //! \param varName Name of variable to be registered
    //! \param ptVal 
    //! \note It is the responsibility of the user to guarantee, that the 
    //!       pointer to the variable is valid throughout the execution. 
    virtual void RegisterExternalVar( unsigned int handle,
                              const std::string& varName,
                              Double * ptVar );
    
    //! Special set function for coordinate values

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
    virtual void SetCoordinates( unsigned int handle,
                         const CoordSystem &coosy,
                         const Vector<Double> &globCoord );
    //@}
    
    // =======================================================================
    //  CALLBACK FUNCTIONALITY
    // =======================================================================
    //@{
        //! \name Callback functionality of mathParser
        
    //! Define signal type of MathParser
    typedef boost::signals2::signal<void()> MathParserSignal;
    
    //! Register callback function for change of value of expression
    boost::signals2::connection 
    virtual AddExpChangeCallBack( const MathParserSignal::slot_function_type
                          &subscriber,
                          unsigned int handle );
    //@}
    
  protected:

    //! Typedef for variable pool (for one handle )
    typedef std::map<std::string, Double>  VarPool;
    
    //! Typedef for mapping handle to variable pool
    typedef std::map<unsigned int, VarPool> PoolMap;
    
    //! Typedef for parser pool
    typedef std::map<unsigned int, mu::Parser> ParserMap;
    
    //! Initialize parser
    void InitParser( mu::Parser &parser, VarPool &pool, 
                     bool isGlobal,
                     bool setDefaults );

    //! Get parser by handle
    mu::Parser & GetParser( unsigned int handle );
    
    //! Factory function for adding variables
    static Double * AddVariable( const char *varName );

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

        
    Double DiffVectorEntry(unsigned int handle, std::string varName, Integer VecPos);

    //! Pool of math parsers
    ParserMap parsers_;
    
    //! Pool of variables
    PoolMap pools_;
    
    //! Map for every global variable parser instance which uses this variable
    std::map<std::string, std::set<unsigned int> > globVarsInUse_;
    
    //! Set with used variables in each expression 
    std::map<unsigned int, std::set<std::string > > varsInUse_;

    //! Set with currently active handles
    std::set<unsigned int> activeHandles_;

    //! Memory for implicitly allocated variables
    static std::list<Double> dynamicPool_;
    
    // =======================================================================
    //  SIGNAL HANDLING STUFF
    // =======================================================================
    
    //! Type definition for shared pointer to MathParserSignal
    typedef shared_ptr<MathParserSignal> PtSig;
    
    //! Callback signals if value of expression changes
    std::map<unsigned int, PtSig > exprChangeSignal_;
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
