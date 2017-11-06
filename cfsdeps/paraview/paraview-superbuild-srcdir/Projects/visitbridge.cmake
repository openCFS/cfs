
#make a dummy visitbridge project so that the user can enable
#or disable building paraview with bridge support
add_external_dummy_project(visitbridge
#  DEPENDS boost
  DEPENDS_OPTIONAL silo cgns
)
