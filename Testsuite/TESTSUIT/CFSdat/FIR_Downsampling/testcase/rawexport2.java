// STAR-CCM+ macro: rawexport2.java
// Written by STAR-CCM+ 11.06.010
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;

public class rawexport2 extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();

    UserFieldFunction userFieldFunction_0 = 
      ((UserFieldFunction) simulation_0.getFieldFunctionManager().getFunction("Variable"));

    userFieldFunction_0.setDefinition("2");

    ImportManager importManager_0 = 
      simulation_0.getImportManager();

    importManager_0.setExportPath("/home/woody/iwaa/iwpa78/cfs/master/Testsuite/TESTSUIT/CFSdat/FIR_Upsampling/testcase/fourcells.case");

    importManager_0.setExportParts(new NeoObjectVector(new Object[] {}));

    Region region_1 = 
      simulation_0.getRegionManager().getRegion("Region");

    Boundary boundary_0 = 
      region_1.getBoundaryManager().getBoundary("Default");

    importManager_0.setExportBoundaries(new NeoObjectVector(new Object[] {boundary_0}));

    importManager_0.setExportRegions(new NeoObjectVector(new Object[] {region_1}));

    importManager_0.setExportScalars(new NeoObjectVector(new Object[] {userFieldFunction_0}));

    importManager_0.setExportVectors(new NeoObjectVector(new Object[] {}));

    importManager_0.setExportOptionAppendToFile(true);

    importManager_0.setExportOptionDataAtVerts(false);

    importManager_0.setExportOptionSolutionOnly(true);

    importManager_0.exportSolution(resolvePath("/home/woody/iwaa/iwpa78/cfs/master/Testsuite/TESTSUIT/CFSdat/FIR_Upsampling/testcase/fourcells.case"), new NeoObjectVector(new Object[] {region_1}), new NeoObjectVector(new Object[] {boundary_0}), new NeoObjectVector(new Object[] {}), new NeoObjectVector(new Object[] {userFieldFunction_0}), true, false);
  }
}
