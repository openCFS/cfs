#ifndef CYL_COORD_SYSTEM_HH
#define CYL_COORD_SYSTEM_HH

#include "coordSystem.hh"

namespace CoupledField {


  //! Describes a cylindrical coordinate system.
  class CylCoordSystem : public CoordSystem{

  public:

    //! Default constructor
    //! \param name (in) name of the coordinate system which is used to 
    //!                  find the additional parameters in the .xml-file
    //! \param ptGrid (in) pointer to finite element grid object
    CylCoordSystem(const std::string & name, Grid * ptGrid);
    
    //! Destructor
    virtual ~CylCoordSystem();

    //! Transform local into global coordinate

    //! This method transforms a point given in local coordinates into a
    //! point w.r.t to the global cartesian x,y,z coordinates.
    //! \param glob (out) point w.r.t. to global cartesian coordinate system
    //! \param loc (in) point w.r.t. to local coordinate system
    void Local2GlobalCoord( Vector<Double> & glob, 
                            const Vector<Double> & loc ) const;
    
    //! Transform global into local coordinate

    //! This method transforms a point given in global cartesian coordinates
    //! into a representation w.r.t to the local coordinate system.
    //! \param loc (out) point w.r.t. to local coordinate system
    //! \param glob (in) point w.r.t. to global cartesian coordinate system
    void Global2LocalCoord( Vector<Double> &loc,  
                            const Vector<Double> & glob ) const;
    
    //! Transform local vector into global one for a given global model point
    
    //! This method transforms a vector with a local coordinate representation
    //! and a given global model point (german: "Aufpunkg") into a vector with 
    //! a representation in global cartesian coordinates.
    //! This method can be used for example to get for a global element 
    //! mid-point and a given local load vector the global representation of
    //! the load vector.
    //! \param globVec (out) vector w.r.t. to global cartesian coordinate 
    //!                      system
    //! \param locVec (in) vector w.r.t. to local coordinate system
    //! \param globModelPoint (in) global model point (german: "Aufpunkt") of
    //!                            of the local/global vector
    void Local2GlobalVector( Vector<Double> & globVec, 
                             const Vector<Double> & locVec, 
                             const Vector<Double> & globModelPoint ) const;
    
    //! Returns for a given coordinate name the according index

    //! This method returns for a given coordinate name (rad,phi,ax)
    //! the according index in the local vector representation.
    //! \param dof (in) name of a coordinate component
    //! \return index of the coordinate component
    UInt GetVecComponent( const std::string & dof ) const;

    //! Returns for a given coordinate index the according name

    //! This method returns for a given coordinate index (1,2,3)
    //! the according name (rad,phi,ax).
    //! \param dof (in) index of the coordinate component
    //! \return name of the coordinate component
    const std::string GetDofName( const UInt dof ) const;


  protected:

    //! Calculate rotation matrix
    void CalcRotationMatrix();
    
    //! Write summary of coordinate system
    void PrintInfo();

    //! second point defining h-axis
    Vector<Double> hAxis_;

    //! second point defining r-axis
    Vector<Double> rAxis_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class CylCoordSystem
  //! 
  //! \purpose 
  //! Thic class implements the representation of clyindrical coordinate
  //! system. It is defined by an point of origin, one point defining the 
  //! the h-axis (3rd component) relative to the origin and one point defining
  //! the r-axis relative to the origin. The latter one is the axis, where 
  //! \f$\rho = 0\f$ holds.
  //! 
  //! \collab 
  //! Ref. to base class CoordSystem.
  //!
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace

#endif
