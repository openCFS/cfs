#ifndef COORD_SYSTEM_HH
#define COORD_SYSTEM_HH

#include <string>

#include "Domain/grid.hh"
#include "Utils/tools.hh"
#include "Utils/vector.hh"
#include "General/environment.hh"

namespace CoupledField {


  //! Base class for describing a coordinate system
  class CoordSystem {

  public:
    
    //! Default constructor
    //! \param name name of the coordinate system which is used to find the
    //!        additional parameters in the .xml-file
    CoordSystem(const std::string & name, Grid * ptGrid);

    //! Destructor
    virtual ~CoordSystem();

    //! Transform local into global coordinate
    virtual void Local2GlobalCoord( Vector<Double> & glob, const Vector<Double> & loc ) = 0;
    
    //! Transform global into local coordinate
    virtual void Global2LocalCoord( Vector<Double> &glob,  const Vector<Double> & glob ) = 0;

    //! Transform a local vector into a global one for a given global point
    virtual void Local2GlobalVector( Vector<Double> & globVec, 
                                     const Vector<Double> & locVec, 
                                     const Vector<Double> & globModelPoint ) = 0;
    

  protected:

    //! helper function for obtaining a point vector
    void GetPoint(Vector<Double> & vec, 
                  StdVector<std::string> & keyVec,
                  StdVector<std::string> & attrVec,
                  StdVector<std::string> & valVec);

    //! Calculate rotation matrix
    virtual void CalcRotationMatrix() = 0;
    
    //! name of the coordinate system
    std::string name_;

    //! Dimension of the coordinate system
    UInt dim_;

    //! Origin of coordinate system
    Vector<Double> origin_;

    //! Rotation matrix (global -> local)
    Matrix<Double> rotationMat_;

    //! Inverse of the rotation matrix (local->global)
    Matrix<Double> invRotationMat_;
    
    //! pointer to grid object
    Grid * ptGrid_;
    
  };


} // end of namespace

#endif
