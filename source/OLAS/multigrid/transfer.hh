// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#ifndef OLAS_TRANSFER_HH
#define OLAS_TRANSFER_HH

//#define DEBUG_TO_CERR
/**********************************************************/
#ifdef DEBUG_TO_CERR
#ifndef DEBUG_TRANSFEROPERATOR
#define DEBUG_TRANSFEROPERATOR
#endif // DEBUG_TRANSFEROPERATOR
#endif // DEBUG_TO_CERR
/**********************************************************/

//#include  "multigrid/totypes.hh"
#include  "multigrid/topology.hh"

namespace CoupledField {
/**********************************************************/

//! encapsulates operators and operations for transfer between AMG levels

/*! This class represents a transfer operator between two AMG levels in
 *  an abstract way. It is created by a call of its member function
 *  CreateOperators. After a successfull call of this function the class
 *  is able to prolongate or interpolate vectors between the AMG levels
 *  or restrict a fine level matrix to the coarse level matrix using the
 *  Galerkin product \f$ A_H = I_h^H A_h I_H^h \f$, with the interpolation
 *  matrix \f$ I_H^h \f$ and its transpose, the restriction matrix
 *  \f$ I_h^H = (I_H^h)^T \f$. <br>
 *  The restriction matrix \f$ I_h^H \f$ is always stored explicitly,
 *  \f$ I_H^h \f$ only, if explicitly demanded (see also CreateOperators).
 *  There exist two implementations for the overloaded function
 *  GalerkinProdukt, that creates the coarse system matrix as Galerkin
 *  product \f$ A_H = I_h^H A_h I_H^h \f$. For the one only the restriction
 *  matrix is needed, for the other the interpolation matrix is used in
 *  addition to the restriction.
 */
  template <typename T>
  class TransferOperator
  {
  public:

    TransferOperator();
    ~TransferOperator();

    //! read only access to prolongation matrix
    const CRS_Matrix<T>* GetProlongation() const {
      return Prolongation_; }
    //! read only access to restriction matrix
    const CRS_Matrix<T>* GetRestriction() const {
      return Restriction_;
    }

    //! resets the class to the state after creation
    void Reset();

    //! creates the prolongation and, if demanded, the restriction

    /*! Creates the restriction matrix and, if demanded the
     *  interpolation matrix. The creation of the interpolation is
     *  not mandatory, but might accelerate the interpolation in return
     *  for occupying again the memory amount of the restriction.
     *  \param matrix the system matrix (so far it is not used,
     *         since we use this very simple interpolation of
     *         unweighted averaging)
     *  \param topology the topology object the C-F splitting
     *         has been created in, it must contain a valid
     *         S graph and a valid splitting
     *  \param itype the interpolation type
     *  \param build_interpolation if true, not only a restriction
     *         matrix \f$ I_H^h \f$ but also a prolongation matrix
     *         \f$ I_H^h := (I_h^H)^T \f$ is built.
     */
    bool CreateOperators( const CRS_Matrix<T>& matrix,
			  const Topology<T>&   topology,
			  const AMGInterpolationType itype,
			  const bool build_interpolation );

    //! constant interpolation
    bool CreateOperatorsConstant( const CRS_Matrix<T>& matrix,
				  const Topology<T>&   topology );
        
    //! simple weighted averaging as interpolation
    bool CreateOperatorsSimpleWeighted( const CRS_Matrix<T>& matrix,
					const Topology<T>&   topology );

    //! constant interpolation, rows of P scaled
    bool CreateOperatorsSmoothedScaling( const CRS_Matrix<T>& matrix,
					 const Topology<T>&   topology );

    //! developing playground for operator building
    bool CreateOperatorsDevelop( const CRS_Matrix<T>& matrix,
				 const Topology<T>&   topology );
        
    //! prolongates the coarse vector v_H to the fine vector v_h

    /*! Prolongates the coarse vector v_H to the fine vector v_h.
     *  Both v_H and v_h must already be set to the matching sizes.
     *  The match of the dimensions (\f$ v_h = I_H^h v_H \f$ with
     *  interpolation \f$ I_H^h \f$) is only checked if the macro
     *  DEBUG_TRANSFEROPERATOR is defined.
     *  \param v_H coarse system vector
     *  \param v_h fine system vector
     *  \param add if true, the result of the prolongation is added
     *         to the vector v_h: \f$ v_h \leftarrow I_H^h v_H \f$
     */
    void Prolongate( const Vector<T>& v_H,
		     Vector<T>& v_h,
		     const bool       add = false ) const;

    //! restricts a fine vector v_h to a coarse vector v_H
    void Restrict( const Vector<T>& v_h,
		   Vector<T>& v_H,
		   const bool       add = false ) const;

    //! creates the coarse system matrix a_H as Galerkin product

    /*! Creates the coarse system matrix a_H as Galerkin product
     *  \f$ A_H = I_h^H A_h I_H^h \f$. This version of the routine
     *  GalerkinProduct does not create temporary matrix data of a_H.
     *  In return for this it needs the graph of \f$ A_H^T \f$
     *  (graph_AHT) and the graph of \f$ (A_h I_H^h)^T \f$ (graph_VT).
     *  These both graphs must be calculated before using the method
     *  Topology<T>::CalcGalerkinGraphs. For the usage of this
     *  method <b>only the restricion</b> matrix must have been built
     *  in this class (see parameter "build_interpolation" of method
     *  CreateOperators). This method expects \f$ A_h \f$ (a_h) to
     *  be <b>sorted conforming to the layout conventions</b>, i.e.
     *  the diagonal entry of each row must be at leading position
     *  followed by the other entries, sorted by column indices.
     *  The coarse matrix built by this routine conforms this layout,
     *  too. See also the second version of this method.
     *  \param a_H the address of a pointer to the coarse matrix.
     *         The method itself will create the matrix object.
     *  \param a_h the matrix of the fine system
     *  \param graph_AHT the graph of \f$ A_H^T \f$
     *  \param graph_VT the graph of \f$ (A_h I_H^h)^T \f$
     */
    void GalerkinProduct(       CRS_Matrix<T>**           a_H,
				const CRS_Matrix<T>&            a_h,
				const DependencyGraph<T>& graph_AHT,
				const DependencyGraph<T>& graph_VT  ) const;

    //! creates the coarse system matrix a_H as Galerkin product

    /*! Creates the coarse system matrix a_H as Galerkin product
     *  \f$ A_H = I_h^H A_h I_H^h \f$. This version of the routine
     *  needs both the interpolation and the restriction matrix to
     *  be created (see parameter "build_interpolation" of method
     *  CreateOperators). In contrary to the other version of this
     *  method this one does not need any matrix graphs precreated,
     *  but will set up the whole coarse matrix in a dynamic class
     *  PreMatrix first and convert it to the CRS format at the end.
     *  This method is the faster alternative! Only, if the creation
     *  of the matrix as PreMatrix in the first step deteriorates
     *  the performance (parallel?), the different implementation
     *  should be chosen. <b>IMPORTANT:</b> This method does not
     *  expect the fine Matrix a_h to conform the layout specifications,
     *  i.e. diagonal entry at leading positions, offdiagonal entries
     *  sorted ascendingly by column indices. It does not require
     *  any of these properties. The created coarse matrix a_H only
     *  fullfills half of these properties, which means that the
     *  diagonal entry can be found at leading position (important
     *  for example for Gauss-Seidel Smoothers), but the offdiagonal
     *  entries are not sorted. (If you want the coarse matrix to be
     *  sorted after creation by this method, change the call of
     *  "preAH.SortDiagonal()" to "preAH.Sort()" in this routine.)
     *  \param a_H address of the pointer, that should keep the
     *         address of the coarse matrix
     *  \param a_h matrix of the fine system
     */
    void GalerkinProduct(       CRS_Matrix<T>** a_H,
				const CRS_Matrix<T>&  a_h ) const;

  protected:

#ifdef TRANSFER_OPERATOR_IMPORT_INTERPOLATION
    void ImportInterpolation( const Integer size_h,
			      const Integer size_H );
#endif

#ifdef DEBUG_TRANSFEROPERATOR
    bool CheckOperators() const;
#endif

    CRS_Matrix<T> *Prolongation_, //!< the interpolation matrix \f$ I_H^h \f$
      *Restriction_;  //!< the restriction matrix \f$ I_h^H \f$
  };

/**********************************************************/
} // namespace CoupledField

#endif // OLAS_TRANSFER_HH
