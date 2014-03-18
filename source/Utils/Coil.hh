#ifndef COIL_FILE_HH
#define COIL_FILE_HH

#include "General/Environment.hh"
#include "MatVec/Vector.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

  // Forward class definition
  class Grid;
  class ParamNode;
  class CoordSystem;


  //! Class for describing coils

  //! This class is basically a container for all parameters required for
  //! describing the different types of coils used in our coupled field
  //! computations.
  //! Besides its container functionality it also provides via its constructor
  //! a method to interface with the XML parameter handling object from which
  //! it obtains the parameter values.
  class Coil {

  public:

    //! Constructor for coils

    //! This is the only allowed constructor for the Coil class. It expects
    //! as first input argument the region identifier of the coil. 
    //! As second input argument
    //! it expects the name of the PDE in whose definition the coil appears.
    //! The constructor will use these two symbolic names to query the XML
    //! parameter handling object for the coil's parameters.
    //! If the coil is a measurement coil it will also, if desired, open the
    //! files for storing the current/voltages and the inductivity.
    Coil( RegionIdType coilRegionId, PtrParamNode coilNode, Grid * ptGrid);

    //! Default desctructor

    //! The default destructor is responsible for closing the output files
    //! in the case that the object describes a measurement coil.
    ~Coil();

    //! Enumeration type for distinguishing the different types of coils
    typedef enum { MEASUREMENT2D, MEASUREMENT3D,
                   VOLTAGE2D, VOLTAGE3D,
                   CURRENT2D, CURRENT3D } Type;

    //! Enumeration type for specifying flow direction for 3D current coils
    typedef enum { NODIR, XDIR, YDIR, ZDIR } FlowDir;

    //! The type of coil this object describes
    Coil::Type coilType_;

    //! The type of coil this object describes as string
    std::string coilTypeS_;

    //! Symbolic name for the coil this object describes
    RegionIdType coilRegionId_;

    //! Name of the region describing the coil
    std::string coilRegionName_;

    //! Area of cross section of a single coil winding
    Double windingCrossSection_;

    //! Size of the voltage (voltage coils) resp. current (current coils)
    std::string value_;

    //! Excitation phase in degrees
    std::string phase_;

    //! Resistance of the coil
    Double resistance_;

    //! Name of results file for inductivity

    //! This string stores the name of the file for writing out the computed
    //! inductivity for the coil. If its value is 'none' this indicates that
    //! the inductivity needs not be computed.
    std::string saveFileL_;

    //! Name of results file for current/voltage

    //! This string stores the name of the file for writing out the computed
    //! voltage for the coil. If its value is 'none' this indicates that the
    //! voltage needs not be computed.
    std::string saveFileU_;

    //! File stream for storing inductivity
    std::ofstream *fileL_;

    //! File stream for storing voltage
    std::ofstream *fileU_;
    
    //! Last step for saving L/U
    Integer lastSaveStep_;

    //! Identifier for windings in 2D coils
    Integer id_;

    //! For current coils the direction of the current flow
    Vector<Double> locFlowDir_;

    //! Reference coordinate system for flow direction
    CoordSystem * flowCoordSys_;

    //! For a current coil in 3D this describes, whether we have a rotational
    //! symmetry alng one axis for the coil.
    bool isRotational_;

    //@{
    //! Component of vector giving midpoint for symmetry axis??
    Double midX_;
    Double midY_;
    Double midZ_;
    //@}

    //@{
    //! Component of vector giving orientation of rotational axis
    Double rotAxisX_;
    Double rotAxisY_;
    Double rotAxisZ_;
    //@}

  private:

    //! The default constructor is not allowed
    Coil();

    //! The copy constructor is not allowed
    Coil( const Coil &c );

  };

}

#endif
