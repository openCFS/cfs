#ifndef POL_COORD_SYSTEM_HH
#define POL_COORD_SYSTEM_HH

#include "CoordSystem.hh"

namespace CoupledField {

  //! Describes a polar coordinate system.
  class PolarCoordSystem : public CoordSystem{

  public:

    //! Default constructor
    //! \param name (in) name of the coordinate system which is used to 
    //!                  find the additional parameters in the .xml-file
    //! \param ptGrid (in) pointer to finite element grid object
    //! \param myParamNode (in) pointer to parameter node of current coosy
    PolarCoordSystem( const std::string & name, Grid * ptGrid, 
                      PtrParamNode myParamNode );
    
    //! Destructor
    virtual ~PolarCoordSystem();
    
    //! Print information about coordinate system to info node
    void ToInfo( PtrParamNode in );
    
    //! Transform local into global coordinate
    
    //! This method transforms a point given in local coordinates into a
    //! point w.r.t to the global cartesian x,y coordinates.
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

    //@{
    //! Transform local vector into global one for a given global model point
    
    //! This method transforms a vector with a local coordinate representation
    //! and a given global model point (german: "Aufpunkt") into a vector with 
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
    
    //! Returns for a given coordinate name the according index

    //! This method returns for a given coordinate name (r,phi)
    //! the according index in the local vector representation.
    //! \param dof (in) name of a coordinate component
    //! \return index of the coordinate component
    UInt GetVecComponent( const std::string & dof, bool silent = false ) const;

    //! Returns for a given coordinate index the according name

    //! This method returns for a given coordinate index (1,2)
    //! the according name (r,phi).
    //! \param dof (in) index of the coordinate component
    //! \return name of the coordinate component
    const std::string GetDofName( const UInt dof ) const;


  protected:

    //! Calculate rotation matrix
    void CalcRotationMatrix();
    
    //! Internal implementation of Local2GlobalVectorInt
    template<class TYPE>
    void Local2GlobalVectorInt( Vector<TYPE> & globVec, 
                                const Vector<TYPE> & locVec, 
                                const Vector<Double> & globModelPoint ) const;
    
    //! second point defining r-axis
    Vector<Double> rAxis_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class PolarCoordSystem
  //! 
  //! \purpose 
  //! Thic class implements the representation of a polar coordinate
  //! system in 2D. It is defined by an point of origin and a vector pointing 
  //! r-direction. The latter one is the axis, where 
  //! \f$\phi = 0\f$ holds.
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
