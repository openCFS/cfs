#ifndef File_NON_CONFORMING_INT
#define File_NON_CONFORMING_INT

#include "baseForm.hh"

namespace CoupledField {

  //! Class for calculation of acoustic-mechanic coupling matrices

  //! This class implements the acoustic-mechanic coupling via surface 
  //! elements for 2D, 2D-axisymmetric and 3D.
  class NonConformingInt : public SurfForm {

  public:
    
    
    //! Constructor
    //! \param dofsPerNode unknowns per node
    //! \param isAxi flag for axisymmetrix setup
    //! \param coordUpdate flag for use of updated lagrangian formulation
    NonConformingInt( UInt dofsPerNode, bool isAxi ,  
                      bool coordUpdate = false  );

    //! Default destructor
    virtual ~NonConformingInt();

    //! Calculate elementwise coupling matrix 

    //! This method calculates the acoustic-mechanic coupling matrix
    //! for a given surface element.
    //! One entry in the coupling matrix \f$C_{pq}\f$ looks like
    //! \f[
    //! C_{pq} = \int_{\Gamma} \rho_0 \left( \begin{array}{c}
    //! \displaystyle N_p N_q n_x \\
    //! \displaystyle N_p N_q n_y \\
    //! \displaystyle N_p N_q n_z \\
    //! \end{array} \right) \; d\Gamma
    //! \f]
    //! where \f$C_{pq}\f$ is a \f$dim \times 1\f$ vector, \f$\rho_0\f$ the
    //! density of water, \f$N_{p,q}\f$ the nodal ansatz functions and 
    //! \f$n_{x,y,z}\f$ the entries of the normal vector \f$\vec{n}\f$
    //! which points outwards from the mechanical domain into the acoustic
    //! one.
    //! The whole element matrix \f$C\f$ is given by
    //! \f[
    //! C = \left( \begin{array}{cccc}
    //! \displaystyle (C_{11}) & (C_{12}) & \cdots & (C_{1n}) \\
    //! \displaystyle (C_{21}) & (C_{22}) & \cdots & (C_{2n}) \\
    //! \displaystyle \cdots & \cdots & \cdots & \cdots \\
    //! \displaystyle (C_{n1}) & (C_{n2}) & \cdots & (C_{nn})
    //! \end{array} \right)
    //! \f]
    //! where \f$n\f$ is the number of nodes of the element. The whole
    //! element matrix \f$C\f$ is of dimension \f$ (n*dim) \times n\f$.
    void CalcElementMatrix( Matrix<Double>& stiffMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 ) ;

  private:
    
    //! Degrees of freedom per node
    UInt dofs_;
  };

} // end of namespace

#endif

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
