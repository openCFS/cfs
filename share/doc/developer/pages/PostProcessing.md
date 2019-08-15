Define a postprocessing result  ([back to main page](/source/CFS_Library_Documentation.md))
===

This section shows which steps are necessary to define a postprocessing result. The steps are illustrated by an example called "Maxwell Force Density". The easiest way is:
* First define your postprocessing result in the xml-Schemas. For our example, the Maxwell Force Density is defined on a surface element. You will find it in [CFS_PDEmagedge.xsd](/share/xml/CFS-Simulation/Schemas/CFS_PDEmagedge.xsd) and it looks like
```
<!-- Definition of surface element result types of magneticedge PDE -->
<xsd:simpleType name="DT_MagneticEdgeSurfElemResult">
  <xsd:restriction base="xsd:token">
    <xsd:enumeration value="magForceMaxwellDensity"/>
  </xsd:restriction>
</xsd:simpleType>
```
* After that, we have to include our result type in our environment enumeration. This you will find in [Environment.cc](/source/General/Environment.cc) as
```
SolutionTypeEnum.Add(MAG_FORCE_MAXWELL_DENSITY, "magForceMaxwellDensity");
```
and add it also in the [Environment.hh](/source/General/Environment.hh) file as
```
// --- flux / derived quantities --
MAG_FLUX_DENSITY, MAG_FLUX, MAG_NORMAL_FLUX_DENSITY, MAG_FLUX_DENSITY_SURF, MAG_FIELD_INTENSITY, MAG_EDDY_CURRENT_DENSITY,
MAG_COIL_CURRENT_DENSITY, MAG_TOTAL_CURRENT_DENSITY, MAG_POTENTIAL_DIV, MAG_FORCE_LORENTZ_DENSITY,
MAG_FORCE_MAXWELL_DENSITY, MAG_EDDY_POWER_DENSITY, MAG_ENERGY_DENSITY, MAG_CORE_LOSS_DENSITY,
MAG_JOULE_LOSS_POWER_DENSITY, MAG_JOULE_POWER_LOSS,
```
* The next step is go into the PDE, section DefinePostProcResults, e.g. in [MagEdgePde.cc](/source/PDE/MagEdgePDE.cc)
* Define a new shared pointer as
```
shared_ptr<ResultInfo> mfd(new ResultInfo);
```
* Allocate some special informations, e.g. the result type, unit of the result, entry type, ...
```
mfd->resultType = MAG_FORCE_MAXWELL_DENSITY;
mfd->dofNames = vecComponents;
mfd->unit = "N/m^3";
mfd->definedOn = ResultInfo::SURF_ELEM;
mfd->entryType = ResultInfo::VECTOR;
```
* The following contains the types of possible results
```
availResults_.insert( mfd );
```
* Now you have to get your results. This can be done for example through a CoefFunction like
```
shared_ptr<CoefFunctionSurfMaxwell> maxForceDens(new CoefFunctionSurfMaxwell(false, matCoefs_, ptGrid_, -1.0, mfd));
```
* Finally, link the ResultInfo and the CoefFunction by
```
DefineFieldResult( maxForceDens, mfd);
```

