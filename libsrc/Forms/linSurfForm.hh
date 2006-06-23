#ifndef FILE_LIN_SURF_FORM
#define FILE_LIN_SURF_FORM

#include "Forms/linearForm.hh"

namespace CoupledField 
{

  //! Base class for all type of right-hand-side surface integrators
  class LinearSurfForm : public LinearForm {

  public:

    //! Standard constructor
    LinearSurfForm();

    //! Destructor
    ~LinearSurfForm();

    //! Set pointer to surface element
    void SetSurfElem( SurfElem * ptSurfElem);

    //! Set normal pointing out of first volume element
    void SetVoluNormal( Vector<Double> & n );
    
    //! Set information of related volume region
    void SetVoluInfo( std::map<RegionIdType, BaseMaterial*> & materials );
    
    //! set additional multiplicative factor for matrix
    void SetFactor(Double factor); 
    
  protected:

    //! Get reference element and coordinates from element iterator
    void ExtractElemInfo( EntityIterator& it);

    //! Current surface element
    const SurfElem * actElem_;

    //! Normal pointing out of first volume element
    Vector<Double> normal_;

    //! Region Ids of associated volume regions
    std::map<RegionIdType, BaseMaterial *> materials_;

    //! Multiplicative factor for vector
    Double factor_;

  };

} // end of namespace
#endif
