#ifndef OLAS_MGPRECOND_HH
#define OLAS_MGPRECOND_HH

/**********************************************************/


#include "OLAS/multigrid/multigrid.hh"

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
class MGPrecond : public BNPrecond< MGPrecond<T>, StdMatrix, T >
{
    public:

        //! entry type of the matrices (e.g. tiny matrices)
        typedef typename AssocType<T>::T_Mtype T_Mtype;
        //! entry type of the vectors (e.g. tiny vectors)
        typedef typename AssocType<T>::T_Vtype T_Vtype;
        //! scalar type (e.g. double, even if T_Mtype is a tiny matrix)
        typedef typename AssocType<T>::T_Stype T_Stype;

        //! Constructor
        //MGPrecond( PtrParamNode params );
        //! second constructor
        MGPrecond( const StdMatrix& matrix,
            PtrParamNode precondNode,
            PtrParamNode olasInfo );

        //! Destructor
        ~MGPrecond();

        //! When called this method returns the type of the preconditioner
        //! object. In the case of an object of this class the return
        //! value is MG (enum constant defined in utils/Environment.hh)
        BasePrecond::PrecondType GetPrecondType() const {
            return BasePrecond::MG;
        };

        //! setup function inherited from class BasePrecond
        void Setup( StdMatrix& sysmatrix );
        
        struct EdgeGeom{
          StdVector<Integer> eNodes; // edge nodes
          Double length;
        };

        //! setup function that takes an additional auxiliary matrix
        void SetupMG(StdMatrix& sysmatrix,
                     StdMatrix& auxmatrix,
                     const AMGType amgType,
                     const StdVector< StdVector< Integer> >& edgeIndNode,
                     const StdVector<Integer>& nodeNumIndex);

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
                    //AMG_->Print( out );
                } else {
                    out << "MGPrecond::Print: AMG preconditioner not"
                           "prepared"<<std::endl;
                }
                return out;
            }

    protected:

  PtrParamNode params_;
  PtrParamNode report_;

  //! the serial solver object
  CoupledField::AMGSolver<T> *AMG_;



    private:

        //! Default constructor

        //! The default constructor is not allowed, since we need size
        //! information and pointers to communication objects for corrected
        //! initialisation.
        MGPrecond() {
            EXCEPTION( "Default constructor of MGPrecond should never be called!");
        };
};

/**********************************************************/
} // namespace CoupledField

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "MGPrecond.cc"
#endif

#endif // OLAS_MGPRECOND_HH
