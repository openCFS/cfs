#include "cylCoordSys.hh"
#include <cmath>

namespace CoupledField{

  CylCoordSystem::CylCoordSystem(const std::string & name,
                                 Grid * ptGrid) 
    : CoordSystem(name, ptGrid) {
    
    ENTER_FCN("CylCoordSystem::CylCoordSystem");
    
    // get origin point of coordinate system
    StdVector<std::string>  keyVec, attrVec, valVec;
    keyVec = "domain", "coordSys", "cylindric", "origin";
    attrVec = "", "", "name";
    valVec = "", "", name;
    GetPoint(origin_, keyVec, attrVec, valVec);

    // get second point for defining the h-axis
    keyVec = "domain", "coordSys", "cylindric", "hAxis";
    attrVec = "", "", "name";
    valVec = "", "", name;
    GetPoint(hAxis_, keyVec, attrVec, valVec);

    // get second point for defining the r-axis
    keyVec = "domain", "coordSys", "cylindric", "rAxis";
    attrVec = "", "", "name";
    valVec = "", "", name;
    GetPoint(rAxis_, keyVec, attrVec, valVec);

    // calculate the rotation matrix and the invers
    CalcRotationMatrix();

    //   // TEST TEST
    //     Vector<Double> globPoint(3), force(3), globVec;
    
    //     force[0] = 0;
    //     force[1] = 1;
    //     force[2] = 45;
    //     globPoint[0] = -4;
    //     globPoint[1] = 8;
    //     globPoint[2] = -1;

    //     Local2GlobalVector(globVec, force, globPoint);

    //     std::cerr << "\n====================\n";
    //     std::cerr << "Global Point is:\n" << globPoint << std::endl;
    //     std::cerr << "Local Force is:\n" << force << std::endl;
    //     std::cerr << "=> Global Force is:\n" << globVec << std::endl;

      
  }
  
  CylCoordSystem::~CylCoordSystem(){
    ENTER_FCN("CylCoordSystem::~CylCoordSystem");
  }

  void CylCoordSystem::Local2GlobalCoord( Vector<Double> & glob, 
                                          const Vector<Double> & loc ) {
    ENTER_FCN("CylCoordSystem::Local2GlobalCoord");
  }
  
  void CylCoordSystem::Global2LocalCoord( Vector<Double> & loc, 
                                          const Vector<Double> & glob ) {
    ENTER_FCN("CylCoordSystem:: Global2LocalCoord");
  }
  
  void CylCoordSystem::
  Local2GlobalVector( Vector<Double> & globVec, 
                      const Vector<Double> & locVec, 
                      const Vector<Double> & globModelPoint ){ 
    ENTER_FCN("CylCoordSystem::Local2GlobalVector");

    Double phi, r;
    static Vector<Double> locModelPoint(3), temp(3), d(3);

    // Transform global cartesian model point into local
    // cartesian one
    rotationMat_.Mult(globModelPoint, locModelPoint);

    std::cerr << "Local Point is " << locModelPoint << std::endl;

    // Calculate the distance and angle to the point
    d =  locModelPoint;
    d -= origin_;

    std::cerr << "d=\n" << d << std::endl;

    r = std::sqrt(d[0] * d[0] + d[1] *d[1]);
    phi = std::atan2(d[1],d[0]);

    std::cerr << "r = " << r << std::endl;
    std::cerr << "phi = " << phi * 180 / PI << std::endl;

    // calculate global coordinate of locVec, by applying the
    // standard conversion routines
    globVec.Resize(3);
    temp.Resize(3);

    // Note: the second term in temp[0] and temp[1] accounts for the
    // tangential component of the locVector
    temp[0] = locVec[0] * std::cos(phi) - locVec[1]*std::sin(phi);
    temp[1] = locVec[0] * std::sin(phi) + locVec[1]*std::cos(phi);
    temp[2] = locVec[2];
    
    std::cerr << "Force before rotation:\n" << temp << std::endl;

    // Now we have the vector in cartesian coordinates for the
    // LOCAL cartesian system. To get the cartesian representation for
    // the GLOBAL one, we have to apply the inverse rotation matrix.
    invRotationMat_.Mult(temp,globVec);

  }
    
  void CylCoordSystem::CalcRotationMatrix() {
    ENTER_FCN( "CalcRotationMatrix" );
        
    Vector<Double> x(3), y(3), z(3);
    Double r_dot;


    // 1) In order to calculate the rotation matrix, we have to find out,
    //    how the local x', y' and 'z' axes are defined, because in a 'global'
    //    cylindric coordinate system we assume, that the z-axis is aligned with
    //    the h-axis and that the axis for rho=0 is equal to the x-axis. Therefore
    //    we can compute the 'local' cartesian axes as follows
    //    z': defined by vector from origin to hAxis point
    //    x': is the normal part of the vector 'origin_-rAxis' w.r.t
    //        z'-Axis
    //    y': implicitly given by cross product of z' and  x' (right hand rule)
    z = hAxis_ - origin_;
    z = z / z.NormL2();

    r_dot = ((rAxis_- origin_) * z) / hAxis_.NormL2();
    x = (rAxis_ - origin_)  - ((hAxis_ *r_dot) / hAxis_.NormL2());
    x = x / x.NormL2();

    y[0] = z[1]*x[2] - z[2]*x[1];
    y[1] = z[2]*x[0] - z[0]*x[2];
    y[2] = z[0]*x[1] - z[1]*x[0];

    // Check if there were coincident points for defining the different axes
    if (x.NormL2() < EPS ||
        y.NormL2() < EPS ||
        z.NormL2() < EPS ) {
      (*error) << "At least two of your points for origin, h-Axis and "
               << "r-Axis coincide in the coordinate system '" << name_
               << "'.\nPlease correct your parameter file!";
      Error( __FILE__, __LINE__ );
    }

    std::cerr << "----------------\n";
    std::cerr << "x_old = \n" << x << std::endl;
    std::cerr << "y_old = \n" << y << std::endl;
    std::cerr << "z_old = \n" << z << std::endl;
    std::cerr << "----------------\n";

    // 2) Calculation of the rotation matrix, which defines the mapping
    //    from the global to the local coordinate system
    //    (ref. Bronstein: Taschenbuch der Mathematik, p. 217f)
    // Note: in order to avoid dividing by zero, an additional check
    //       is performed, if the x/y/z-component is in the order of
    //       machine precission.

    rotationMat_.Resize(3,3);
    rotationMat_[0][0] = (std::abs(x[0]) < EPS ) ? 0.0 : (1.0 / x[0]);
    rotationMat_[0][1] = (std::abs(x[1]) < EPS ) ? 0.0 : (1.0 / x[1]);
    rotationMat_[0][2] = (std::abs(x[2]) < EPS ) ? 0.0 : (1.0 / x[2]);

    rotationMat_[1][0] = (std::abs(y[0]) < EPS ) ? 0.0 : (1.0 / y[0]);
    rotationMat_[1][1] = (std::abs(y[1]) < EPS ) ? 0.0 : (1.0 / y[1]);
    rotationMat_[1][2] = (std::abs(y[2]) < EPS ) ? 0.0 : (1.0 / y[2]);

    rotationMat_[2][0] = (std::abs(z[0]) < EPS ) ? 0.0 : (1.0 / z[0]);
    rotationMat_[2][1] = (std::abs(z[1]) < EPS ) ? 0.0 : (1.0 / z[1]);
    rotationMat_[2][2] = (std::abs(z[2]) < EPS ) ? 0.0 : (1.0 / z[2]);

    // 3) Calculate invers rotation matrix, which defines mapping from
    //    local to global cartesian coordinate system
    rotationMat_.Invert(invRotationMat_);

    std::cerr << "Rotation Matrix:\n" << rotationMat_ << std::endl;

    
  }

} // end of namespace
