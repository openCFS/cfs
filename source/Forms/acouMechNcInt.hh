#ifndef _ACOU_MECH_NC_INT_2007_
#define _ACOU_MECH_NC_INT_2007_

#include "General/defs.hh"
#include "General/environment.hh"
#include "baseForm.hh"

namespace CoupledField {

  //! Class for calculation of non-conforming grid acoustic-mechanic coupling
  //! matrices

  //! This class implements the acoustic-mechanic coupling via non-conforming 
  //! grid interfaces for 2D, 2D-axissymmetric and 3D.
class EntityIterator;
template <class TYPE> class Matrix;

  class AcouMechNcInt : public SurfForm {

  public:
    
    
    //! Constructor
    //! \param dofsPerNode unknowns per node
    //! \param isAxi flag for axisymmetric setup
    AcouMechNcInt( UInt dofsPerNode, bool isAxi, bool coordUpdate = false );

    //! Default destructor
    virtual ~AcouMechNcInt();

    // Set formulation of integrator (ACOU_POTENTIAL/ACOU_PRESSURE)
    void SetFormulation(SolutionType aformulation);

    //! Calculate elementwise coupling matrix 

    /*! This method calculates the acoustic-mechanic coupling matrix
        for a given intersection(!) element.
        One entry in the coupling matrix \f$C_{pq}\f$ looks like
        \f[
        C_{pq} = \int_{\Gamma} \rho_0 \left( \begin{array}{c}
        \displaystyle N_p N_q n_x \\
        \displaystyle N_p N_q n_y \\
        \displaystyle N_p N_q n_z \\
        \end{array} \right) \; d\Gamma
        \f]
        where \f$C_{pq}\f$ is a \f$dim \times 1\f$ vector, \f$\rho_0\f$ the
        density of the fluid and \f$n_{x,y,z}\f$ the entries of the normal
        vector \f$\vec{n}\f$ which points outwards from the mechanical domain
        into the acoustic one. \f$N_{p,q}\f$ are the nodal ansatz functions of
        the two parent elements on master and slave side of the interface.
        The whole element matrix \f$C\f$ is given by
        \f[
        C = \left( \begin{array}{cccc}
        \displaystyle (C_{11}) & (C_{12}) & \cdots & (C_{1n}) \\
        \displaystyle (C_{21}) & (C_{22}) & \cdots & (C_{2n}) \\
        \displaystyle \cdots & \cdots & \cdots & \cdots \\
        \displaystyle (C_{n1}) & (C_{n2}) & \cdots & (C_{nn})
        \end{array} \right)
        \f]
        where \f$n\f$ is the number of nodes of the element. The whole
        element matrix \f$C\f$ is of dimension \f$ (n*dim) \times n\f$.
    */
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );

  private:
    
    //! Degrees of freedom per node
    UInt dofs_;

    //! Formulation of acouPDE (acouPotential | acouPressure)
    SolutionType formulation_;

  };
} // namespace CoupledField

#endif // _ACOU_MECH_NC_INT_2007_
