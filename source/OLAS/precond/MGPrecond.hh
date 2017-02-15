#ifndef OLAS_MGPRECOND_HH
#define OLAS_MGPRECOND_HH

/**********************************************************/

#include "multigrid/multiGrid.hh"

#include "BasePrecond.hh"
#include "BNPrecond.hh"

namespace CoupledField {
/**********************************************************/

//! Base Class for Multigrid Preconditioners

/*! Base class for Multigrid Preconditioners. <b>For CRS-Matrices only
 *  </b>. This class is basically a wrapper for class AMGSolver. For
 *  a detailed documentations of the AMG preconditioner and the
 *  corresponding settings in OLAS_Params see the documentation of class
 *  AMG_Solver.
 */
template <typename T>
class MGPrecond :
public BNPrecond< MGPrecond<T>, StdMatrix, T >
{
    public:
        //! Constructor
        MGPrecond( OLAS_Params* params = NULL );
        //! second constructor
        MGPrecond( const StdMatrix&   matrix,
                   OLAS_Params* params,
                   OLAS_Report* report );
        //! Destructor
        ~MGPrecond();

        //! When called this method returns the type of the preconditioner
        //! object. In the case of an object of this class the return
        //! value is MG (enum constant defined in utils/Environment.hh)
        PrecondType GetPrecondType() const {
            return MG;
        };
        
        //! setup function inherited from class BasePrecond
        void Setup( StdMatrix& sysmatrix );
        
        //! setup function that takes an additional auxiliary matrix
        void Setup( const StdMatrix& sysmatrix,
                    const StdMatrix& auxmatrix );

        //! preconditioning step, inherited from class BasePrecond

        /*! Performs one AMG preconditioning step. The parameter
         *  \c sysmatrix is useless in this derivation, since the
         *  address of the matrix, that is used, has been saved during
         *  the setup process.
         *  \param sysmatrix \b ignored
         *  \param rhs right hand side of the linear equation
         *  \param sol solution (\b output only, will be initialized first)
         */
        void Apply( const StdMatrix& sysmatrix,
                    const SingleVector& rhs,
                          SingleVector& sol ) const;

        //! print method for the AMG preconditioner
		std::ostream& Print( std::ostream& out ) const {
            if( AMG_ ) {
                AMG_->Print( out );
            } else {
                out << "MGPrecond::Print: AMG preconditioner not"
                       "prepared"<<std::endl;
            }    
            return out;
        }

    protected:

		OLAS_Params *params_;
		OLAS_Report *report_;
    	//! the serial solver object
        AMGSolver<T>    *AMG_;

    private:

        //! Default constructor

        //! The default constructor is not allowed, since we need size
        //! information and pointers to communication objects for corrected
        //! initialisation.
        MGPrecond() {
            Error( "Default constructor of MGPrecond should never be called!",
                   __FILE__, __LINE__ );
        };
};

/**********************************************************/
} // namespace CoupledField

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "mgprecond.cc"
#endif

#endif // OLAS_MGPRECOND_HH
