#ifndef CYL_COORD_SYSTEM_HH
#define CYL_COORD_SYSTEM_HH

#include "coordSystem.hh"

namespace CoupledField {


  //! Describes a cylindrical coordinate system.
  //! \note At the moment the the cylindrical system has to be aligned
  //! with the global x,y, and z axis, since no rotation of the 
  //! coordinate system is performed.
  class CylCoordSystem : public CoordSystem{

  public:

    //! Default constructor
    //! \param name name of the coordinate system which is used to find the
    //!        additional parameters in the .xml-file
    CylCoordSystem(const std::string & name, Grid * ptGrid);
    
    //! Destructor
    virtual ~CylCoordSystem();

    //! Transform local into global coordinate
    void Local2GlobalCoord( Vector<Double> & glob, const Vector<Double> & loc );
    
    //! Transform global into local coordinate
    void Global2LocalCoord( Vector<Double> & loc, const Vector<Double> & glob );

    //! Transform a local vector into a global one for a given global point
    void Local2GlobalVector( Vector<Double> & globVec, 
                             const Vector<Double> & locVec, 
                             const Vector<Double> & globModelPoint );
    

  protected:

    //! Calculate rotation matrix
    void CalcRotationMatrix();
    
    //! second point defining h-axis
    Vector<Double> hAxis_;

    //! second point defining r-axis
    Vector<Double> rAxis_;

  };

} // end of namespace

#endif
