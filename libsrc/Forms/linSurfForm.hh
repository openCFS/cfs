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
    void SetVoluInfo( const StdVector<RegionIdType> & regionIds,
                      const StdVector<BaseMaterial*>& materials );
    
    //! set additional multiplicative factor for matrix
    void SetFactor(Double factor); 
    
  protected:

    //! Current surface element
    SurfElem * actElem_;

    //! Normal pointing out of first volume element
    Vector<Double> normal_;

    //! Region Ids of associated volume regions
    StdVector<RegionIdType> regionIds_;

    //! Materials of volume regions
    StdVector<BaseMaterial*> materials_;

    //! Multiplicative factor for vector
    Double factor_;

  };

} // end of namespace
#endif
