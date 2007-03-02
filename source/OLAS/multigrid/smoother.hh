/* $Id$ */

#ifndef OLAS_SMOOTHER_HH
#define OLAS_SMOOTHER_HH

/**********************************************************/

#include "matvec/matvec.hh"
#include "multigrid/topology.hh"

namespace OLAS {
/**********************************************************/

//! base class for multigrid smoothers

/*! The base class for all smoothers in a multigrid hierarchy.
 *
 */
template <typename T>
class Smoother
{
    public:

        //! entry type of the matrices (e.g. tiny matrices)
        typedef typename assocType<T>::T_Mtype T_Mtype;
        //! entry type of the vectors (e.g. tiny vectors)
        typedef typename assocType<T>::T_Vtype T_Vtype;
        //! scalar type (e.g. double, even if T_Mtype is a tiny matrix)
        typedef typename assocType<T>::T_Stype T_Stype;

        //! enumeration type for the smoothing direction
        enum Direction { FORWARD = 0, BACKWARD };

        //! constructor
        Smoother();
        //! destructor
        virtual ~Smoother();

        //! simply "return prepared_;"
        bool IsPrepared() const { return prepared_; }

        //! eventual setup of the smoother

        /*! This function prepares the smoother, if necessary.
         *  In derived classes this method must set prepared_ = true.
         *  \param matrix the problem matrix
         */
        virtual bool Setup( const CRS_Matrix<T>& matrix );

        //! One smoothing step.

        /*! One smoothing step of the problem Ax = b.
         *  \param matrix The problem matrix A
         *  \param rhs Vector with the right hand side b
         *  \param sol Vector, where the solution x should be stored
         *         into
         *  \param direction If a non-symmetric smoother is applied in
         *         the same way for pre-smoothing and post-smoothing,
         *         the AMG preconditioner is far from being symmetric.
         *         This asymmetry can be reduced by using "mirrored"
         *         versions of a smoother for pre- and post-smoothing
         *         respectively. For example in a Gauss-Seidel smoother
         *         this can be achieved by updating the unknowns once
         *         in a ascending order for pre-smoothing and in
         *         descending order for postsmoothing.\n
         *         This update direction is triggered by the parameter
         *         \c direction. Since there might be different ways of
         *         using the "mirrored" variants than just using the one
         *         version for pre-smoothing the other for post-smoothing
         *         (e.g. alternating application in case of more than
         *         one smoothing step), the constant name is not assotiated
         *         to the application as PRESMOOTHER or POSTSMOOTHER,
         *         just in case you wondered about the names. The actual
         *         implementations in derived classes might be different
         *         from just using different orderings, or e.g. also
         *         ignore this parameter, like damped Jacobi.
         *  \param force_setup If the smoother object is not prepared
         *         (see IsPrepared(), Setup()), Setup() is called before
         *         the smoothing step. If the setup has already been done,
         *         Setup() ist not called, except force_setup is set to
         *         true. To make sure that the smoother works correctly,
         *         the default value of force_setup is true, for savety
         *         reasons. Set force_setup to true, if it is granted that
         *         the setup can be reused.
         */
        virtual void Step( const CRS_Matrix<T>& matrix,
                           const Vector<T>&     rhs,
                                 Vector<T>&     sol,
                           const Direction      direction   = FORWARD,
                           const bool           force_setup = true ) = 0;


        //! resets the object to the state after creation

        /*! Resets the smoother object to the state after creation.
         *  Call Smoother::Reset() in the Reset() of derived classes.
         */
        virtual void Reset();

    protected:

        //! constant reference to the parameter object
        
        //! flag that signs, whether the smoother has been prepared
        bool prepared_;
};

/**********************************************************/
} // namespace OLAS

#endif // OLAS_SMOOTHER_HH
