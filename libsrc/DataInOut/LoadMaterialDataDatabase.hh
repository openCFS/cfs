#ifndef FILE_LOADMATERIALDATADATABASE
#define FILE_LOADMATERIALDATADATABASE

/**************************************************************************/
/* File:   LoadMaterialDataDatabase.hh                                            */
/* Author: Fred Hofer                                                     */
/* Date:   14. Dez. 99                                                    */
/*                                                                        */
/* Reads Material Data from a Materialfile.                               */
/**************************************************************************/


#include <Utils/vector.hh>
#include <General/environment.hh>
#include "DataInOut/Database/database.hh"
#include "MaterialData.hh"

namespace CoupledField
{

//! Class for reading material data file
/*! 
  Includes all functions for reading material data out of file
*/
class LoadMaterialDataDatabase
{
private:
  //! Establish connection to database
  void InitializeDatabase();

  //! Read piezo material data
  void ReadPiezo (MaterialData &matData, const std::string &matName);


  //! Read piezo material linear stiffness data
  void ReadStiffnessLinear (MaterialData &matData, const std::string &matName, int matidx);

  //! Read piezo material nonlinear stiffness data f(S) - not yet implemented
  void ReadStiffnessNonlinear_S(MaterialData &matData, const std::string &matName, int matidx){Error("Nonlinear piezo matrix not yet implemented",__FILE__,__LINE__);};

  //! Read piezo material stiffness data, hysteresis f(S)  - not yet implemented
  void ReadStiffnessHysteresis_S(MaterialData &matData, const std::string &matName, int matidx){Error("Piezo hysteresis not yet implemented",__FILE__,__LINE__);};


  //! Read piezo material linear coupling data
  void ReadCouplingLinear (MaterialData &matData, const std::string &matName, int matidx);

  //! Read piezo material nonlinear coupling data f(E) - not yet implemented
  void ReadCouplingNonlinear_E(MaterialData &matData, const std::string &matName, int matidx){Error("Nonlinear coupling not yet implemented",__FILE__,__LINE__);};

  //! Read piezo material nonlinear coupling data f(S) - not yet implemented
  void ReadCouplingNonlinear_S(MaterialData &matData, const std::string &matName, int matidx){Error("Nonlinear coupling not yet implemented",__FILE__,__LINE__);};

  //! Read piezo material coupling data, hysteresis f(S)  - not yet implemented
  void ReadCouplingHysteresis_E(MaterialData &matData, const std::string &matName, int matidx){Error("Coupling hysteresis not yet implemented",__FILE__,__LINE__);};

  //! Read piezo material coupling data, hysteresis f(S)  - not yet implemented
  void ReadCouplingHysteresis_S(MaterialData &matData, const std::string &matName, int matidx){Error("Coupling hysteresis not yet implemented",__FILE__,__LINE__);};


  //! Read piezo material linear dielectric data
  void ReadDielectricLinear (MaterialData &matData, const std::string &matName, int matidx);

  //! Read piezo material nonlinear coupling data f(E) - not yet implemented
  void ReadDielectricNonlinear_E(MaterialData &matData, const std::string &matName, int matidx){Error("Nonlinear dielectricity not yet implemented",__FILE__,__LINE__);};

  //! Read piezo material coupling data, hysteresis f(S)  - not yet implemented
  void ReadDielectricHysteresis_E(MaterialData &matData, const std::string &matName, int matidx){Error("Dielectricity hysteresis not yet implemented",__FILE__,__LINE__);};


  //! Read material linear conductivity data
  void ReadConductivityLinear (MaterialData &matData, const std::string &matName, int matidx);

  //! Read material nonlinear conductivity data f(T) - not yet implemented
  void ReadConductivityNonlinear_T(MaterialData &matData, const std::string &matName, int matidx){Error("Nonlinear conductivity not yet implemented",__FILE__,__LINE__);};

  //! Read material nonlinear conductivity data f(B) - not yet implemented
  void ReadConductivityNonlinear_B(MaterialData &matData, const std::string &matName, int matidx){Error("Nonlinear conducitvity not yet implemented",__FILE__,__LINE__);};

  //! Read material conductivity data, hysteresis f(T)  - not yet implemented
  void ReadConductivityHysteresis_T(MaterialData &matData, const std::string &matName, int matidx){Error("Conductivity hysteresis not yet implemented",__FILE__,__LINE__);};

  //! Read material conductivity data, hysteresis f(B)  - not yet implemented
  void ReadConductivityHysteresis_B(MaterialData &matData, const std::string &matName, int matidx){Error("Conductivity hysteresis not yet implemented",__FILE__,__LINE__);};


  //! Read material linear permeability data
  void ReadPermeabilityLinear (MaterialData &matData, const std::string &matName, int matidx);

  //! Read material nonlinear permeability data f(T) - not yet implemented
  void ReadPermeabilityNonlinear_T(MaterialData &matData, const std::string &matName, int matidx){Error("Nonlinear permeability not yet implemented",__FILE__,__LINE__);};

  //! Read material nonlinear permeability data f(B) - not yet implemented
  void ReadPermeabilityNonlinear_B(MaterialData &matData, const std::string &matName, int matidx){Error("Nonlinear permeability not yet implemented",__FILE__,__LINE__);};

  //! Read material permeability data, hysteresis f(T)  - not yet implemented
  void ReadPermeabilityHysteresis_T(MaterialData &matData, const std::string &matName, int matidx){Error("Permeability hysteresis not yet implemented",__FILE__,__LINE__);};

  //! Read material permeability data, hysteresis f(B)  - not yet implemented
  void ReadPermeabilityHysteresis_B(MaterialData &matData, const std::string &matName, int matidx){Error("Permeability hysteresis not yet implemented",__FILE__,__LINE__);};


  //! Read material linear dielectric data
  void ReadMagnetisationLinear (MaterialData &matData, const std::string &matName, int matidx);

  //! Read material nonlinear magnetisation data f(B) - not yet implemented
  void ReadMagnetisationNonlinear_B(MaterialData &matData, const std::string &matName, int matidx){Error("Nonlinear magnetisation not yet implemented",__FILE__,__LINE__);};

  //! Read material magnetisation data, hysteresis f(B)  - not yet implemented
  void ReadMagnetisationHysteresis_B(MaterialData &matData, const std::string &matName, int matidx){Error("Magnetisation hysteresis not yet implemented",__FILE__,__LINE__);};



  //! Read fluid material data
  void ReadFluid (MaterialData &matData, const std::string &matName);

  //! Read magnetic material data
  void ReadMagnetic (MaterialData &matData, const std::string &matName);

  //! Database connection
  Database Db_;

public:
  //! Default constructor
  LoadMaterialDataDatabase();

  //! Overloaded constructor for the same use as in LoadMaterial and LoadMaterialFile
  LoadMaterialDataDatabase(const char * fileName);

  //! Default destructor (to close database connection)
  ~LoadMaterialDataDatabase();

  /// private constructor, loads all materials defined in the parserfile from the materialfile "filename"
  void GetMaterial( MaterialData& material, const std::string matName, const std::string matType);
  
};

} // end namespace CoupledField
 
#endif
