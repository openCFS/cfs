// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id: AFWsmoother.hh 15666 2017-05-12 14:41:17Z kroppert $ */

#ifndef OLAS_AFWSMOOTHER_HH
#define OLAS_AFWSMOOTHER_HH

/**********************************************************/

#include "OLAS/multigrid/smoother.hh"

namespace CoupledField {
/**********************************************************/

//! Arnold-Falk-Winther AFW-smoother for edge-AMG

/*! AFW-smoother for edge-AMG.
 *  This class implements a block-wise edge smoother, where
 *  all edges, connected by one common node are smoothed together
 */

template <typename T>
class AFWSmoother : public Smoother<T>
{
    public:

        AFWSmoother();
        ~AFWSmoother();

        //@{
        //! get value for the SOR damping parameter Omega
        Double& GetOmega() { return Omega_; }
        const Double& GetOmega() const { return Omega_; }

        //! set value for the SOR damping parameter Omega
        void SetOmega( const Double& omega ) {
            Omega_ = omega;
        }
        //@}
        

        //! Create the connection of node-edges which are smoothed together

        /*! We store for each node i in the auxiliary matrix a set
         * of srcMatrix-indices (edges) which are connected to node i.
         */
        void CreatePatches(const CRS_Matrix<Double>& AuxMatrix,
                           const StdVector< StdVector< Integer> >& edgeIndNode,
                           const StdVector<Integer> nodeNumIndex);


        //! setup of the AFW smoother

        /*! This function prepares the AFW smoother.
         *  That means that an array is created, containing the
         *  inverses of the diagonal entries of the matrix "matrix".
         *  (reimplementation of the virtual function of the base
         *  class Smoother
         *  \param matrix system matrix
         *  \return true, if the setup was successful
         */
        bool Setup( const CRS_Matrix<T>& matrix );

        //! one Gauss-Seidel/SOR step

        /*! One Gauss-Seidel/SOR step
         *  \param matrix The problem matrix A
         *  \param rhs Vector with the right hand side b
         *  \param sol Vector, where the solution x should be stored
         *         into
         *  \param direction \c Smoother::FORWARD: the unknowns are
         *         updated in ascending order starting at 1.
         *         \c Smoother::BACKWARD: vice versa. For a general
         *         description of this parameter see Smoother::Step.
         *  \param force_setup flag for the control of the setup
         *         (in class AFWSmoother the creation of array containing
         *         inverse values of the matrix' diagonal entries, for
         *         further information of this parameter see
         *         Smoother::Setup)
         */
        void Step( const CRS_Matrix<T>& matrix,
                   const Vector<T>&     rhs,
                         Vector<T>&     sol,
                   const typename Smoother<T>::Direction direction
                                        = Smoother<T>::FORWARD,
                   const bool           force_setup = true );
        
        //! resets the object to the state after its creation
        void Reset();


    private:
        //! Create the extraction-blocks of the system matrix

        /*! Extract blocks from the system matrix and rhs-vector (at least
         * store the needed indices) and perform a local inversion of this
         * extracted matrix (corresponds to the node-wise inversion of the diagonal
         * element in the classical Gauss-Seidel.
         * Every extracted matrix corresponds to a patch as defined
         * in Patches_. This means every node is assigned such a small inverted
         * matrix and a set of indices in the rhs-vector.
         * This is just the preparation for the Setup-method
         */
        void ExtractPatches( const CRS_Matrix<T>& SysMat);

    protected:

        UInt Size_; //!< size of the system
        UInt SizeNodes_;  //!< number of nodes in the auxiliary-matrix
        Double Omega_; //!< damping factor for a SOR

        StdVector< StdVector< UInt> > Patches_;

        //! contains the inverted extraction matrix for every node,
        //! as defined in Patches_
        StdVector< Matrix<T> > InvExtMat_;

        //! contains the extracted matrix for every node,
        //! as defined in Patches_
        StdVector< Matrix<T> > ExtMat_;

        //! contains the indices which have to be extracted from
        //! the rhs-vector
        StdVector< StdVector<UInt> > ExtRhs_;

    private:
        const UInt *pRow_;
        const UInt *pCol_;
        const T *pDat_;
};

/**********************************************************/
} // namespace CoupledField

#endif // OLAS_AFWSMOOTHER_HH
