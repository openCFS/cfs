#include "DataInOut/LoadMaterialDataDatabase.hh"
#include "General/environment.hh"
#include "WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{

LoadMaterialDataDatabase::LoadMaterialDataDatabase()
{
  ENTER_FCN("LoadMaterialDataDatabase::LoadMaterialDataDatabase()");
  InitializeDatabase();
}

LoadMaterialDataDatabase::LoadMaterialDataDatabase(const char* fileName)
{
  ENTER_FCN("LoadMaterialDataDatabase::LoadMaterialDataDatabase(const char*)");
  InitializeDatabase();
}

void LoadMaterialDataDatabase::InitializeDatabase()
{
  ENTER_FCN("LoadMaterialDataDatabase::InitializeDatabase");

  std::string hostName, userName, passwd, matDatabaseName;
  Integer port;
  
  // get connection parameters
  params->Get("hostName", hostName, "database");
  params->Get("port", port, "database");
  params->Get("materialDatabaseName", matDatabaseName, "database");
  params->Get("userName", userName, "database");
  params->Get("password", passwd, "database");
  
  Db_.Connect (hostName, port, userName,
	       passwd, matDatabaseName);
}

LoadMaterialDataDatabase::~LoadMaterialDataDatabase()
{
  ENTER_FCN("LoadMaterialDataDatabase::~LoadMaterialDataDatabase");
  Db_.Close();
}

void LoadMaterialDataDatabase::GetMaterial(MaterialData &material, 
                                           const std::string matName, 
                                           const std::string matType)
{
  ENTER_FCN("LoadMaterialDataDatabase::GetMaterial");
#ifdef DEBUG
  (*debug)<<"Loading material data of "<<matName<<" from type "<<matType<<" out of database"<<std::endl;
#endif
  if (matType=="piezo")
    ReadPiezo(material,matName);
  else if (matType=="fluid")
    ReadFluid(material,matName);
  else if (matType=="magnetic")
    ReadMagnetic(material,matName);
  else
    Error("Material not found in database");
} 

void LoadMaterialDataDatabase::ReadPiezo (MaterialData &matData, const std::string &matName)
{
  ENTER_FCN("LoadMaterialDataDatabase::ReadPiezo");
  matData.DefFull3dMatrix();
  // Get line with general material data
  std::string where("mat_name='");
  where += matName; where+="'";
  Db_.SelectFrom("idx,mat_name,mechanic,coupling,dielectric,density,alpha,beta","Material_parameter",where);
  dbMatrix matParam;
  if (!Db_.FetchFields(matParam))
    Error("LoadMaterialDataDatabase::ReadPiezo : FetchFields(matParam)");
  if (matParam.getNoOfRow()>1)
    Warning("More than one entry found for material");
  else if (matParam.getNoOfRow()<1)
  {
    std::string err1("Material "); err1+=matName; err1+=" not found in database.";
    Error(err1.c_str());
  }
  int matidx;
  matParam["idx"]->get(matidx,0);

  // Mechanic stiffness
  std::string mechdecider;
  matParam["mechanic"]->get(mechdecider,0);
  if (mechdecider=="uninitialized")
    Error ("Try to access uninitialized material parameter 'mechanic'");
  else if (mechdecider=="linear")
    ReadStiffnessLinear(matData,matName,matidx);
  else if (mechdecider=="nonlinear_S")
    ReadStiffnessNonlinear_S(matData,matName,matidx);
  else if (mechdecider=="hysteresis_S")
    ReadStiffnessHysteresis_S(matData,matName,matidx);
  else
    Error ("Unknown result at field 'mechanic' in table 'Material_parameter'");

  // Coupling
  std::string coupdecider;
  matParam["coupling"]->get(coupdecider,0);
  if (coupdecider=="uninitialized")
    Error ("Try to access uninitialized material parameter 'coupling'");
  else if (coupdecider=="linear")
    ReadCouplingLinear(matData,matName,matidx);
  else if (coupdecider=="nonlinear_E")
    ReadCouplingNonlinear_E(matData,matName,matidx);
  else if (coupdecider=="nonlinear_S")
    ReadCouplingNonlinear_S(matData,matName,matidx);
  else if (coupdecider=="hysteresis_E")
    ReadCouplingHysteresis_E(matData,matName,matidx);
  else if (coupdecider=="hysteresis_S")
    ReadCouplingHysteresis_S(matData,matName,matidx);
  else 
    Error ("Unknown result at field 'coupling' in table 'Material_parameter'");
  
  // Dielectric
  std::string dieldecider;
  matParam["dielectric"]->get(dieldecider,0);
  if (dieldecider=="uninitialized")
    Error ("Try to access uninitialized material parameter 'dielectric'");
  else if (dieldecider=="linear")
    ReadDielectricLinear(matData,matName,matidx);
  else if (dieldecider=="nonlinear_E")
    ReadDielectricNonlinear_E(matData,matName,matidx);
  else if (dieldecider=="hysteresis_E")
    ReadDielectricHysteresis_E(matData,matName,matidx);
  else 
    Error ("Unknown result at field 'dielectric' in table 'Material_parameter'");

  Double lambda, mu;
  matData.GetPiezoMatrixData(0,0,mu);
  matData.GetPiezoMatrixData(0,1,lambda);
  matData.SetLambda(lambda);
  matData.SetMu(mu);

  Double density, alpha, beta;
  matParam["density"]->get(density,0);
  matParam["alpha"]->get(alpha,0);
  matParam["beta"]->get(beta,0);
  matData.SetDensity(density);
  matData.SetDampingCoeffs(alpha,beta);

  matData.SetName(matName.c_str());

  Info->PrintPiezoMat(matData);
}

void LoadMaterialDataDatabase::ReadStiffnessLinear(MaterialData &matData, 
                                                   const std::string &matName, 
                                                   int matidx)
{
  ENTER_FCN("LoadMaterialDataDatabase::ReadStiffnessLinear");
  std::stringstream wheremech;
  wheremech<<"material_idx="<<matidx;
  std::string fields("v11,v12,v13,v14,v15,v16,v21,v22,v23,v24,v25,v26,v31,v32,v33,v34,v35,v36,v41,v42,v43,v44,v45,v46,");
  fields += "v51,v52,v53,v54,v55,v56,v61,v62,v63,v64,v65,v66";
  Db_.SelectFrom(fields,"Material_mech_lin",wheremech.str());
  dbMatrix mechmatrix;
  if (!Db_.FetchFields(mechmatrix))
    Error("LoadMaterialDataDatabase::ReadStiffnessLinear : FetchFields(mechmatrix)");
  if (mechmatrix.getNoOfRow()>1)
    Warning("More than one entry found for mechanical stiffness matrix for material");
  else if (mechmatrix.getNoOfRow()<1)
    Error ("No mechanical stiffness entry found for material");
  std::stringstream fieldindex;
  double val;
  dbColumn *col;
  int i,j;
  for (i=1; i<=6; i++)
  {
    for (j=1; j<=6; j++)
    {
      fieldindex.str("");
      fieldindex<<"v"<<i<<j;
      col = mechmatrix[fieldindex.str()];//fieldindex.str()];
      col->get(val,0);
      matData.SetPiezoMatrixData(i-1,j-1,val);
    } 
  }
}

void LoadMaterialDataDatabase::ReadCouplingLinear (MaterialData &matData, 
                                                   const std::string &matName, 
                                                   int matidx)
{
  ENTER_FCN("LoadMaterialDataDatabase::ReadCouplingLinear");
  // Get piezoelectric coupling terms
  std::stringstream wherecoupling;
  wherecoupling<<"material_idx="<<matidx;
  std::string fields;
  fields = "v11,v12,v13,v14,v15,v16,v21,v22,v23,v24,v25,v26,v31,v32,v33,v34,v35,v36";
  Db_.SelectFrom(fields,"Material_coupling_lin",wherecoupling.str());
  dbMatrix couplingmatrix;
  if (!Db_.FetchFields(couplingmatrix))
    Error("LoadMaterialDataDatabase::ReadCouplingLinear FetchFields(couplingmatrix");
  if (couplingmatrix.getNoOfRow()>1)
    Warning("More than one entry found for mechanical stiffness matrix for material");
  else if (couplingmatrix.getNoOfRow()<1)
    Error("No piezoelectric coupling term found for material");
  std::stringstream fieldindex;
  dbColumn *col;
  double val;
  for (int i=1; i<=3; i++)
  {
    for (int j=1; j<=6; j++)
    {
      fieldindex.str("");
      fieldindex<<"v"<<i<<j;
      col = couplingmatrix[fieldindex.str()];
      col->get(val,0);
      matData.SetPiezoMatrixData(i+6-1,j-1,val);
      matData.SetPiezoMatrixData(j-1,i+6-1,val);
    }
  }
}

void LoadMaterialDataDatabase::ReadDielectricLinear (MaterialData &matData, 
                                                     const std::string &matName, 
                                                     int matidx)
{
  ENTER_FCN("LoadMaterialDataDatabase::ReadDielectricLinear");
  // Get dielectric terms
  std::stringstream wheredielec;
  wheredielec<<"material_idx="<<matidx;
  std::string fields;
  fields = "v11,v12,v13,v21,v22,v23,v31,v32,v33";
  Db_.SelectFrom(fields,"Material_dielec_lin",wheredielec.str());
  dbMatrix dielecmatrix;
  if (!Db_.FetchFields(dielecmatrix))
    Error("LoadMaterialDataDatabase::ReadDielectricLinear FetchFields(dielecmatrix");
  if (dielecmatrix.getNoOfRow()>1)
    Warning("More than one entry found for mechanical stiffness matrix for material");
  else if (dielecmatrix.getNoOfRow()<1)
    Error("No dielectric term found for material");
  std::stringstream fieldindex;
  dbColumn *col;
  double val;
  for (int i=1; i<=3; i++)
  {
    for (int j=1; j<=3; j++)
    {
      fieldindex.str("");
      fieldindex<<"v"<<i<<j;
      col = dielecmatrix[fieldindex.str()];
      col->get(val,0);
      matData.SetPiezoMatrixData(i+6-1,j+6-1,val);
    }
  }
}


void LoadMaterialDataDatabase::ReadFluid (MaterialData &matData, const std::string &matName)
{
  ENTER_FCN("LoadMaterialDataDatabase::ReadFluid");
  // Get line with general material data
  std::string where("mat_name='");
  where += matName; where+="'";
  Db_.SelectFrom("idx,mat_name,density,compression,alpha,beta","Material_parameter",where);
  dbMatrix matParam;
  if (!Db_.FetchFields(matParam))
    Error("LoadMaterialDataDatabase::ReadFluid : FetchFields(matParam)");
  if (matParam.getNoOfRow()>1)
    Warning("More than one entry found for material");
  else if (matParam.getNoOfRow()<1)
  {
    std::string err1("Material "); err1+=matName; err1+=" not found in database.";
    Error(err1.c_str());
  }
  
  double compress, density,alpha,beta;
  dbColumn *col;
  col = matParam["compression"];
  col->get(compress,0);
  matData.SetCompressibility(compress);
  col = matParam["density"];
  col->get(density,0);
  matData.SetDensity(density);
  col = matParam["alpha"];
  col->get(alpha,0);
  col = matParam["beta"];
  col->get(beta,0);
  matData.SetDampingCoeffs(alpha,beta);
  
  matData.SetName(matName.c_str());
  Info->PrintFluidMat(matData);
}


void LoadMaterialDataDatabase::ReadMagnetic (MaterialData &matData, const std::string &matName)
{
  ENTER_FCN("LoadMaterialDataDatabase::ReadMagnetic");
  std::string where("mat_name='");
  where += matName; where+="'";
  Db_.SelectFrom("idx,conductivity,permeability,magnetisation","Material_parameter",where);
  dbMatrix matParam;
  if (!Db_.FetchFields(matParam))
    Error("LoadMaterialDataDatabase::ReadMagnetic : FetchFields(matParam)");
  if (matParam.getNoOfRow()>1)
    Warning("More than one entry found for material");
  else if (matParam.getNoOfRow()<1)
  {
    std::string err1("Material "); err1+=matName; err1+=" not found in database.";
    Error(err1.c_str());
  }
  int matidx;
  matParam["idx"]->get(matidx,0);
  
  // Conductivity
  std::string conducdecider;
  matParam["conductivity"]->get(conducdecider,0);
  if (conducdecider=="uninitialized")
    Error ("Try to access uninitialized material parameter 'conductivity'");
  else if (conducdecider=="linear")
    ReadConductivityLinear(matData,matName,matidx);
  else if (conducdecider=="nonlinear_T")
    ReadConductivityNonlinear_T(matData,matName,matidx);
  else if (conducdecider=="nonlinear_B")
    ReadConductivityNonlinear_B(matData,matName,matidx);
  else if (conducdecider=="hysteresis_T")
    ReadConductivityHysteresis_T(matData,matName,matidx);
  else if (conducdecider=="hysteresis_B")
    ReadConductivityHysteresis_B(matData,matName,matidx);
  else
    Error ("Unknown result at field 'conductivity' in table 'Material_parameter'");

  // Permeability
  std::string permeadecider;
  matParam["permeability"]->get(permeadecider,0);
  if (permeadecider=="uninitialized")
    Error ("Try to access uninitialized material parameter 'permeability'");
  else if (permeadecider=="linear")
    ReadPermeabilityLinear(matData,matName,matidx);
  else if (permeadecider=="nonlinear_T")
    ReadPermeabilityNonlinear_T(matData,matName,matidx);
  else if (permeadecider=="nonlinear_B")
    ReadPermeabilityNonlinear_B(matData,matName,matidx);
  else if (permeadecider=="hysteresis_T")
    ReadPermeabilityHysteresis_T(matData,matName,matidx);
  else if (permeadecider=="hysteresis_B")
    ReadPermeabilityHysteresis_B(matData,matName,matidx);
  else
    Error ("Unknown result at field 'permeability' in table 'Material_parameter'");

  // Magnetisation
  std::string magnetdecider;
  matParam["magnetisation"]->get(magnetdecider,0);
  if (magnetdecider=="uninitialized")
    Error ("Try to access uninitialized material parameter 'magnetisation'");
  else if (magnetdecider=="linear")
    ReadMagnetisationLinear(matData,matName,matidx);
  else if (magnetdecider=="nonlinear_B")
    ReadMagnetisationNonlinear_B(matData,matName,matidx);
  else if (magnetdecider=="hysteresis_B")
    ReadMagnetisationHysteresis_B(matData,matName,matidx);
  else
    Error ("Unknown result at field 'magnetisation' in table 'Material_parameter'");
  matData.SetName(matName.c_str());
  Info->PrintMagMat(matData);
}

void LoadMaterialDataDatabase::ReadConductivityLinear(MaterialData &matData, 
                                                      const std::string &matName, 
                                                      int matidx)
{
  ENTER_FCN("LoadMaterialDataDatabase::ReadConductivityLinear");
  std::stringstream whereconduc;
  whereconduc<<"material_idx="<<matidx;
  std::string fields;
  fields = "v11,v12,v13,v21,v22,v23,v31,v32,v33";
  Db_.SelectFrom(fields,"Material_conduc_lin",whereconduc.str());
  dbMatrix conducmatrix;
  if (!Db_.FetchFields(conducmatrix))
    Error("LoadMaterialDataDatabase::ReadConductivityLinear FetchFields(conducmatrix");
  if (conducmatrix.getNoOfRow()>1)
    Warning("More than one entry found for conductivity matrix for material");
  else if (conducmatrix.getNoOfRow()<1)
    Error("No conductivity term found for material");
  std::stringstream fieldindex;
  dbColumn *col;
  double val;
  for (int i=1; i<=3; i++)
  {
    for (int j=1; j<=3; j++)
    {
      fieldindex.str("");
      fieldindex<<"v"<<i<<j;
      col = conducmatrix[fieldindex.str()];
      col->get(val,0);
      matData.SetConductivity(i-1,j-1,val);
    }
  }
}

void LoadMaterialDataDatabase::ReadPermeabilityLinear(MaterialData &matData, 
                                                      const std::string &matName, 
                                                      int matidx)
{
  ENTER_FCN("LoadMaterialDataDatabase::ReadPermeabilityLinear");
  std::stringstream whereperm;
  whereperm<<"material_idx="<<matidx;
  std::string fields;
  fields = "v11,v12,v13,v21,v22,v23,v31,v32,v33";
  Db_.SelectFrom(fields,"Material_perm_lin",whereperm.str());
  dbMatrix permeamatrix;
  if (!Db_.FetchFields(permeamatrix))
    Error("LoadMaterialDataDatabase::ReadPermeabilityLinear FetchFields(permeamatrix");
  if (permeamatrix.getNoOfRow()>1)
    Warning("More than one entry found for permeability matrix for material");
  else if (permeamatrix.getNoOfRow()<1)
    Error("No permeability term found for material");
  std::stringstream fieldindex;
  dbColumn *col;
  double val;
  for (int i=1; i<=3; i++)
  {
    for (int j=1; j<=3; j++)
    {
      fieldindex.str("");
      fieldindex<<"v"<<i<<j;
      col = permeamatrix[fieldindex.str()];
      col->get(val,0);
      matData.SetPermeability(i-1,j-1,val);
    }
  }
}

void LoadMaterialDataDatabase::ReadMagnetisationLinear (MaterialData &matData, 
                                                        const std::string &matName, 
                                                        int matidx)
{
  ENTER_FCN("LoadMaterialDataDatabase::ReadMagnetisationLinear");
  std::stringstream where;
  where<<"material_idx="<<matidx;
  std::string fields;
  fields = "v1,v2,v3";
  Db_.SelectFrom(fields,"Material_magnet_lin",where.str());
  dbMatrix magnetmatrix;

  if (!Db_.FetchFields(magnetmatrix))
    Error("LoadMaterialDataDatabase::ReadMagnetisationLinear FetchFields(permeamatrix");
  if (magnetmatrix.getNoOfRow()>1)
    Warning("More than one entry found for magnetisation matrix for material");
  else if (magnetmatrix.getNoOfRow()<1)
    Error("No magnetisation term found for material");
  std::stringstream fieldindex;
  dbColumn *col;
  double mX,mY,mZ;
  fieldindex.str("v1");
  col = magnetmatrix[fieldindex.str()];
  col->get(mX,0);
  fieldindex.str("v2");
  col = magnetmatrix[fieldindex.str()];
  col->get(mY,0);
  fieldindex.str("v3");
  col = magnetmatrix[fieldindex.str()];
  col->get(mZ,0);
  matData.SetPermMag(mX,mY,mZ);
}

} // End of namespace
