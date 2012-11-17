MACRO(ADD_OPTION
    OPT_NAME_LONGFORM
    OPT_NAME_DATATYPE
    OPT_NAME_DEFAULT
    OPT_NAME_SHORTHELP
    ADDITIONAL_LATEX_HELP
    )

  STRING(REGEX REPLACE "SEMICOLON" ";" SHORTHELP_CC ${OPT_NAME_SHORTHELP})

  # Generate C++ source
  IF("${OPT_NAME_DATATYPE}" STREQUAL "string")
    STRING(CONFIGURE "#define CPLREADER_OPT_${OPT_NAME_LONGFORM} CPLR_POM(\"${OPT_NAME_LONGFORM}\", \\
                                                   std::${OPT_NAME_DATATYPE}, \\
                                                   param_${OPT_NAME_LONGFORM}, \\
                                                   \"${OPT_NAME_DEFAULT}\", \\
                                                   \"${SHORTHELP_CC}\")" OPT_DEFINE)
  ELSE("${OPT_NAME_DATATYPE}" STREQUAL "string")
    STRING(CONFIGURE "#define CPLREADER_OPT_${OPT_NAME_LONGFORM} CPLR_POM(\"${OPT_NAME_LONGFORM}\", \\
                                                   ${OPT_NAME_DATATYPE}, \\
                                                   param_${OPT_NAME_LONGFORM}, \\
                                                   ${OPT_NAME_DEFAULT}, \\
                                                   \"${SHORTHELP_CC}\")" OPT_DEFINE)
  ENDIF("${OPT_NAME_DATATYPE}" STREQUAL "string")

  SET(OPT_DEFINES "${OPT_DEFINES}\n${OPT_DEFINE}")

  IF("${OPT_NAME_DATATYPE}" STREQUAL "string")
    STRING(CONFIGURE
      "std::${OPT_NAME_DATATYPE} param_${OPT_NAME_LONGFORM} = \"${OPT_NAME_DEFAULT}\";"
      VARIABLE_NAME)
  ELSE("${OPT_NAME_DATATYPE}" STREQUAL "string")
    STRING(CONFIGURE
      "${OPT_NAME_DATATYPE} param_${OPT_NAME_LONGFORM} = ${OPT_NAME_DEFAULT};"
      VARIABLE_NAME)
  ENDIF("${OPT_NAME_DATATYPE}" STREQUAL "string")

  SET(VARIABLE_NAMES "${VARIABLE_NAMES}    ${VARIABLE_NAME}\n")

  IF("${OPT_NAME_DATATYPE}" STREQUAL "string")
    SET(SETTING "settings.SetString(\"${OPT_NAME_LONGFORM}\", param_${OPT_NAME_LONGFORM});")
  ELSE("${OPT_NAME_DATATYPE}" STREQUAL "string")
    IF("${OPT_NAME_DATATYPE}" STREQUAL "double")
      STRING(CONFIGURE "settings.SetDouble(\"${OPT_NAME_LONGFORM}\", param_${OPT_NAME_LONGFORM});" SETTING)
    ELSE("${OPT_NAME_DATATYPE}" STREQUAL "double")
      STRING(CONFIGURE "settings.SetInt(\"${OPT_NAME_LONGFORM}\", param_${OPT_NAME_LONGFORM});" SETTING)
    ENDIF("${OPT_NAME_DATATYPE}" STREQUAL "double")
  ENDIF("${OPT_NAME_DATATYPE}" STREQUAL "string")

  SET(SETTINGS "${SETTINGS}    ${SETTING}\n")

  SET(BOOST_OPTIONS "${BOOST_OPTIONS}        CPLREADER_OPT_${OPT_NAME_LONGFORM}\n")

  # Generate entries for LaTex longtable
  STRING(REGEX REPLACE "_" "\\\\_" SHORTHELP ${OPT_NAME_SHORTHELP})
  STRING(REGEX REPLACE "\\|" "\\\\textbar\\\\," SHORTHELP ${SHORTHELP})
  STRING(REGEX REPLACE "SPACE" "\\\\textvisiblespace\\\\," SHORTHELP ${SHORTHELP})
  STRING(REGEX REPLACE "SEMICOLON" ";" SHORTHELP ${SHORTHELP})
  SET(LATEX_OPTIONS "${LATEX_OPTIONS}    
        -\\,-${OPT_NAME_LONGFORM} & ${OPT_NAME_DEFAULT} & ${SHORTHELP}\\\\    
        & & \\multicolumn{1}{p{120mm}}{    
          ${ADDITIONAL_LATEX_HELP}
        }\\\\    
        \\hline\n"
    )

ENDMACRO(ADD_OPTION)


ADD_OPTION(docu
  string
  ""
  "Create directory with detailed documentation (can just be PDF)"
  "This option will create the directory {\\\\tt cplreader\\\\_docs} in the
   current working directory and put the file {\\\\tt cplreader\\\\_docu.pdf} there."
  )

ADD_OPTION(name
  string
  "" 
  "Name of dataset as well as of the directory containing the dataset"
  "This is the name of the dataset to be worked on. For most readers
   the data files must be placed in a directory of this name, too."
  )

ADD_OPTION(type
  string
  CFX
  "Type of dataset (can be ANSYS | CFS++ | CFX | CFX_EXPORT | FASTEST | OPENFOAM | CGNS | ENSIGHT | FIELDVIEW)"
  "Specify the type of dataset to be read: ANSYS see Sec. \\\\ref{sec:ansys},
   CFX see Sec. \\\\ref{sec:cfx}, CFX\\\\_EXPORT see Sec. \\\\ref{sec:cfx-export},
   FASTEST see Sec. \\\\ref{sec:fastest}, OPENFOAM see Sec. \\\\ref{sec:openfoam}."
  )

ADD_OPTION(xmlfile
  string
  ""
  "XML file containing options for cplreader."
  "The XML file specified here contains informations for the file readers/writers."
  )

ADD_OPTION(basedir
  string
  .
  "Base directory where subdirectories are searched."
  "The files to be read are searched for in the directory specified by basedir.

   For CFX this means for example:\\\\newline\\\\newline
   {\\\\tt
   basedir/name/\\\\newline
   \\\\textbar-- name\\\\newline
   \\\\textbar \\\\quad  \\\\textbar-- 1.trn\\\\newline
   \\\\textbar \\\\quad  +-- 2.trn\\\\newline
   \\\\textbar-- name.def\\\\newline
   +-- name.res\\\\newline
   }
   \\\\newline
   By default the files will be searched for in the current working directory."
  )

ADD_OPTION(extfiles
  bool
  true
  "Use external time step files (just like CFX .trn files)."
  "If this parameter is true \\\\cplreader\\\\, will create on master .h5 file containing
   the mesh and just references to the time step files. For each time step an extra
   {\\\\it external} .h5 file will be created. If this parameter is false everything
   goes into one big .h5 file. It may be tricky to copy such big files to USB hard
   disks, which are usually formatted with a FAT32 file system."
  )
  
ADD_OPTION(pressureRhsForWave
  bool
  false
  "Compute wave equation source term based on Laplacian of fluid pressure."
  "If this parameter is true, the acouRhsLoad for the wave equation is computed by the Lighthill
   Tensor but by the Laplacian of the fluid pressure."
  )

ADD_OPTION(dim
  uint32_t
  3
  "Dimension of fluid data. (can be 2 | 3)"
  "This parameter defines the geometrical dimension of the fluid data set.
   Most readers do not require this parameter. It is however present for
   the FASTEST reader."
  )

ADD_OPTION(numsteps
  uint32_t
  0
  "Number N of timesteps to read. Default: read all available steps."
  "The number of timesteps to read. By default (0) all available timesteps will be read."
  )

ADD_OPTION(activeparts
  string
  all
  "Values will only be output on partitions specified by activeparts (all | numbers seperated by SPACE or SEMICOLON or |)."
  "When dealing with very large data sets with many partitions/sub-domains, it
   is often convenient to read the fields and compute the acoustic source for
   just a few of them."
  )

ADD_OPTION(outputfields
  string
  acouRhsLoad
  "Which physical fields should be output ([acouDivLighthillTensor | acouRhsLoad | acouRhsLoadDensity | fluidMechDensity | fluidMechPressure | fluidMechVelocity | fluidMechTKE | fluidMechSkinFriction | acouLambVec | acouLambRhs | aeroAcouSourceRhs]).  Values may be separated by SPACE or SEMICOLON or |"
  "Sometimes it may be necessary to do some post-processing on the fluid fields
   to get some understanding of the problem at hand. For this reason it is possible
   to write the most important fields (velocity, pressure, source terms, turbulence
   kinetic energy, etc.) to the .h5 files by using the outputfields parameter."
  )

ADD_OPTION(firststep
  uint32_t
  1
  "Index of first time step (only for CFX_EXPORT, FASTEST)"
  "To start reading at a different time step than the first one this parameter can by used.
   It is just implemented for CFX\\\\_EXPORT and FASTEST at the moment."
  )

ADD_OPTION(stepincr
  uint32_t
  1
  "Step increment for reading the files"
  "Increment time step counter by this value und just read each stepincr-th time step."
  )

ADD_OPTION(timestep
  double
  -1.0
  "Time step length T in seconds"
  "This parameter specifies the length of each timestep in seconds."
  )
  
ADD_OPTION(calcMeanPresField
  bool
  0
  "Toggle calculation of mean pressure field"
  "This parameter toggles the calculation of the time averaged mean pressure field"
  )
  
ADD_OPTION(calcMeanVelField
  bool
  0
  "Toggle calculation of mean velocity field"
  "This parameter toggles the calculation of the time averaged mean velocity field"
  )
  
    
ADD_OPTION(numStepsForMeanField
  uint32_t
  0
  "Specify how many files to consider for the mean field calculation"
  "By default we will take the value given by the numSteps parameter to perform the 
   temporal average of the flow and pressure field. By this parameter we can take more or 
   less steps for the mean flow computation."
  )

ADD_OPTION(calcsrc
  bool
  false
  "Calculate the acoustic sources from velocity"
  "If this parameter is set a velocity field has to present in each time step from which
   the Lighthill source term will be computed."
  )

ADD_OPTION(geomscale
  string
  "1.0 1.0 1.0"
  "Scale factors for geometry in the format factorX factorY factorZ"
  "This parameter specifies the factors to be multiplied with the corresponding
   components of the nodal coordinates. This can be useful, if the CFD computation
   used a different scaling (cf. OpenFOAM simulations of the Freiberg guys)"
  )

ADD_OPTION(velscale
  string
  "1.0 1.0 1.0"
  "Scale factors for velocity field in the format factorX factorY factorZ"
  "This parameter specifies the factors to be multiplied with the corresponding
   components of the velocity field. This can be useful, if the CFD computation
   used a different scaling (cf. OpenFOAM simulations of the Freiberg guys)"
  )
  
ADD_OPTION(presscale
  string
  "1.0"
  "Scale factors for pressure field in the format factor"
  "This parameter specifies the factor to be multiplied with the pressure field. 
   This can be useful, if the CFD computation
   used a different scaling (cf. OpenFOAM simulations of the Freiberg guys)"
  )

ADD_OPTION(density
  double
  1.0
  "Density of the fluid. Default density is 1.0 for testing."
  "The mean density of the fluid may be specified by this parameter.
   Some practically relevant values include:
   air $1.2041 kg/m^3$ (at $20^\\\\circ C$ and $101.325 kPa$), water $998.2 kg/m^3$ (at $20^\\\\circ C$)"
  )

ADD_OPTION(segfault
  bool
  true
  "Cause program to segfault if an exception is encountered."
  "This parameter is especially important for debugging purposes and
   should not be necessary in standard runs."
  )

ADD_OPTION(deltemp
  bool
  true
  "Delete temporary files."
  "Some file readers (ANSYS) create or use temporary files. This parameter determines,
   whether those files will be deleted or not."
  )

ADD_OPTION(outprec
  string
  double
  "Write data as single or as double precision (can be single | double)"
  "With outprec the user can specify the precision with which the floating point numbers
   will be written to the .h5 files. specifying {\\\\it single} here may save some space
   but the user has to be sure what he is doing!"
  )

ADD_OPTION(justinit
  bool
  false
  "Just perform initialization step."
  "Specifying justinit will let \\\\cplreader\\\\, just perform the init phase.
   This comes in handy when developing a new reader."
  )

ADD_OPTION(justmesh
  bool
  false
  "Just convert the mesh."
  "Specifying justmesh will make \\\\cplreader\\\\, just read the mesh without any other data."
  )

ADD_OPTION(degen
  bool
  false
  "Allow degenerated element shapes."
  "This option just applies to the ANSYS reader. If true all elements other than lines will
   degenerated version of quadrilaterals and hexahedrons."
  )

ADD_OPTION(strict
  bool
  false
  "Be strict."
  "For the ANSYS reader this option means that \\\\cplreader\\\\, will complain about holes
   in the node list and connectivity, i.e. if the user forgot to nummrg the entities in ANSYS.
   For the CFX reader this option means that if corrupt transient files have been found, the
   reader will complain and stop working."
  )

ADD_OPTION(verbose
  bool
  false
  "Be verbose"
  "If this parameter is set to true lots of helpful information will be output which is
   not necessery to have in most cases and would normally just clutter the shell."
  )

ADD_OPTION(deffile
  string
  ""
  "Definition file name. Only for CFX!"
  "Non-standard path to .def file for CFX."
  )

ADD_OPTION(trntol
  double
  0.05
  "Ratio by which .trn files may differ."
  "Ratio  by   which  .trn   files  may  differ   before  they   are  regarded
  corrupt.  E.g.:  Let's  assume the  mean  size  of  .trn  files is  1MB.  By
  specifying --trntol  0.05 the files may  vary in size from  0.95MB to 1.05MB
  before an error condition  is assumed."
  )

ADD_OPTION(lhsrc
  string
  ""
  "Column of Lighthill source term in FASTEST result files (e.g. col1)."
  "If FASTEST already calculated the Lighthill source term, the column in the ASCII file
   can be specified using this parameter."
  )

ADD_OPTION(vx
  string
  ""
  "Column of x-velocity in FASTEST result files (e.g. col2)."
  "Specify column of x-velocity in FASTEST ASCII files."
  )

ADD_OPTION(vy
  string
  ""
  "Column of y-velocity in FASTEST result files (e.g. col3)."
  "Specify column of y-velocity in FASTEST ASCII files."
  )

ADD_OPTION(vz
  string
  ""
  "Column of z-velocity in FASTEST result files (e.g. col4)."
  "Specify column of z-velocity in FASTEST ASCII files."
  )

ADD_OPTION(pres
  string
  ""
  "Column of pressure in FASTEST result files (e.g. col5)."
  "Specify column of pressure in FASTEST ASCII files."
  )

ADD_OPTION(reduceElementOrder
  bool
  false
  "disregard high order nodes"
  "If this parameter is set, higher order noes of the elements will be
  disregarded for the calclation. Only reasonable for quadratic elements."
  )

ADD_OPTION(doIntAverageCentre
  bool
  false
  "averages the integration over nodes at centre"
  "If this parameter is set acouRhsLoad is calculated by using the velocity
  vector at the centre of each element. May help with quadratic elements."
  )

ADD_OPTION(doCalcMultiNodes
  bool
  true
  "This flag needs to be set to remedy region interface artifacts."
  "If multiple regions exist, nodes which are shared by both are stored multiple
  times, once for each region. But their acouRhsLoad is only calculated in each
  region seperatly and therefore get different values which need to be added.
  This flag should be set on otherwise artifacts at the region interface occur."
  )

ADD_OPTION(calcSurfaceTerms
  string
  none
  "This flag enables the calculation of the surface integral within for lighthill source term. Specify 'all' or give a list of surface regions separated by SPACE or SEMICOLON or |"
  "This flag indicates if we wish to calculate the surface integral of the 
  Lighthill right hand side which is up to now a very costly operation and therefore 
  off by default."
  )

ADD_OPTION(cfxUseStnFrame
  bool
  false
  "Use velocity in stationary frame (CFX only)"
  "When moving meshes are used in CFX, the velocity is usually computed in
  the moving frame (i.e. velocity relative to the moving mesh). Activate this
  option if you need the velocity in the stationary frame (i.e. absolute
  velocity)."
  )
 
ADD_OPTION(cfxSinglePrecision
  string
  "auto"
  "Single or double precision datasets. 'auto' for auto-detection. (CFX only)"
  "Tells cplreader if it should expect float or double precision datasets in CFX result files. Set to 'auto' for auto-detection."
  )

ADD_OPTION(cfxLastStep
  int32_t
  -1
  "Last time step computed by CFX simulation run"
  "Specifies the number of the last time step computed by the CFX simulation
  run (not necessarily the same as the last step that was written out). This
  option must be given for data produced by CFX 12 and later. For data
  produced by CFX 11 and earlier this option can be omitted, because it is
  determined automatically."
  )
   
CONFIGURE_FILE("ParamsInit.cc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/ParamsInit.cc")

CONFIGURE_FILE("Documentation/latex/cplreader_options.tex.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Documentation/latex/cplreader_options.tex")
