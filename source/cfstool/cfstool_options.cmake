MACRO(ADD_OPTION
    OPT_NAME_LONGFORM
    OPT_NAME_SHORTFORM
    OPT_NAME_DATATYPE
    OPT_NAME_DEFAULT
    OPT_NAME_SHORTHELP
    ADDITIONAL_LATEX_HELP
    )

  STRING(REGEX REPLACE "SEMICOLON" ";" SHORTHELP_CC ${OPT_NAME_SHORTHELP})
  IF(NOT ${OPT_NAME_SHORTFORM} STREQUAL "")
    SET(OPT_NAME  "${OPT_NAME_LONGFORM},${OPT_NAME_SHORTFORM}")
  ELSE(NOT ${OPT_NAME_SHORTFORM} STREQUAL "")
    SET(OPT_NAME  "${OPT_NAME_LONGFORM}")
  ENDIF(NOT ${OPT_NAME_SHORTFORM} STREQUAL "")


  # Generate C++ source
  IF("${OPT_NAME_DATATYPE}" STREQUAL "string")
    STRING(CONFIGURE "#define CFSTOOL_OPT_${OPT_NAME_LONGFORM} \\
                      CFSTOOL_POM(\"${OPT_NAME}\", \\
                                                   std::${OPT_NAME_DATATYPE}, \\
                                                   param_${OPT_NAME_LONGFORM}, \\
                                                   \"${OPT_NAME_DEFAULT}\", \\
                                                   \"${SHORTHELP_CC}\")" OPT_DEFINE)
  ELSE("${OPT_NAME_DATATYPE}" STREQUAL "string")
    STRING(CONFIGURE "#define CFSTOOL_OPT_${OPT_NAME_LONGFORM} \\
                      CFSTOOL_POM(\"${OPT_NAME}\", \\
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
    SET(PARAMNODE "leafs.Push_back(PtrParamNode(new ParamNode()));
    leafs[leafs.GetSize()-1]->SetName(\"${OPT_NAME_LONGFORM}\");
    leafs[leafs.GetSize()-1]->SetValue(param_${OPT_NAME_LONGFORM});")
  ELSE("${OPT_NAME_DATATYPE}" STREQUAL "string")
    SET(PARAMNODE "leafs.Push_back(PtrParamNode(new ParamNode()));
    leafs[leafs.GetSize()-1]->SetName(\"${OPT_NAME_LONGFORM}\");
    leafs[leafs.GetSize()-1]->SetValue(lexical_cast<std::string>(param_${OPT_NAME_LONGFORM}));")
  ENDIF("${OPT_NAME_DATATYPE}" STREQUAL "string")

  SET(PARAMNODES "${PARAMNODES}    ${PARAMNODE}\n")

  SET(BOOST_OPTIONS "${BOOST_OPTIONS}        CFSTOOL_OPT_${OPT_NAME_LONGFORM}\n")

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

ADD_OPTION(files f
  string
  ""
  "choose reference_file compare_file and if required outfile"
  ""
  )

ADD_OPTION(eps e
  string
  1.0e-6
  "Sets the tolerance for check if two vector are the same."
  ""
  )

ADD_OPTION(freeCoord ""
  string
  "" 
  "Arguments needed: 'comp start stop inc'.\\\\n\\\\
E.g.  'x 0 1 0.1'"
  ""
  )

ADD_OPTION(mode m
  string
  scalardiff
  "Type of mode:  \\\\n\\\\
calcAverage | convert | scalardiff |meshdiff | meshdiffnormed"
  ""
  )

ADD_OPTION(forceSegFault ""
  string
  "false"
  "force a segmentation fault at exceptions"
  ""
  )


CONFIGURE_FILE("ParamsInit.cc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/ParamsInit.cc")

# TODO: Documentation
#CONFIGURE_FILE("Documentation/latex/cplreader_options.tex.in"
#  "${CMAKE_CURRENT_BINARY_DIR}/Documentation/latex/cplreader_options.tex")
