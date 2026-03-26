#ifndef DEFAULT_COORD_SYSTEM_HH
#define DEFAULT_COORD_SYSTEM_HH

#include "CoordSystem.hh"

namespace CoupledField {


  //! Describes the global cartesian coordinate system
  class DefaultCoordSystem : public CoordSystem{

  public:

    //! Default constructor
    //! \param name (in) name of the coordinate system which is used to 
    //!                  find the additional parameters in the .xml-file
    DefaultCoordSystem(Grid * ptGrid );
    
    //! Destructor
    virtual ~DefaultCoordSystem();

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

    //@{
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
    void Local2GlobalVector( Vector<Complex> & globVec, 
                             const Vector<Complex> & locVec, 
                             const Vector<Double> & globModelPoint ) const;
    //@}
    

    //! Returns the global rotation matrix for a given point

    //! This method returns the rotation matrix defining defining a rotation,
    //! by which the global coordinate system has to be rotated, so that
    //! so that it represents the current one in that point.
    //! \param rotMatrix rotation matrix for global point
    //! \param point point w.r.t. to global Cartesian coordinate system
    virtual void 
    GetGlobRotationMatrix( Matrix<Double> & rotMatrix,
                           const Vector<Double>& point ) const;

    //! Returns the full 3x3 global rotation matrix for a given point

    //! This method returns the full 3x3 rotation matrix defining defining 
    //! a rotation, by which the global coordinate system has to be rotated, 
    //! it represents the current one in that point.
    //! \param rotMatrix rotation matrix for global point
    //! \param point point w.r.t. to global Cartesian coordinate system
    virtual void 
    GetFullGlobRotationMatrix( Matrix<Double> & rotMatrix,
                               const Vector<Double>& point ) const;

    /** For 2D "x"->1 and "y"->2 but in the axisymmetric case also "r"->1 and "z"->2
     *  For 3D "z"->3 
     * @see CoordSystem::GetVecComponent() */
    UInt GetVecComponent( const std::string & dof, bool silent = false ) const;

    /** 1-> "x", 2-> "y" and for 3D 3-> "z"
     * aysymmetric case is not considered but could be implemented */
     const std::string GetDofName( const UInt dof ) const;


  protected:

    //! Calculate rotation matrix
    void CalcRotationMatrix(){};
    
  };


#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class DefaultCoordSystem
  //! 
  //! \purpose 
  //! This class mainly a dummy class, since it represents the global
  //! cartesian coordinate system, so no real transformation is applied.
  //! Its main purpose is to convert dofs and the according indices.
  //! 
  //! \collab 
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
