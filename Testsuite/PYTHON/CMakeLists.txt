# this section of test-cases tests Python pre- and postprocessing scriptsg
if(TESTSUITE_PYTHON)
  # test via doctest
  add_subdirectory(doctest/hdf5_tools)
  add_subdirectory(doctest/nrfcommon)
  add_subdirectory(doctest/nrfreader)
  add_subdirectory(doctest/tutorial)
  # this section of test-cases tests the visualization routines implemented in matviz.py
  add_subdirectory(matviz/tensor)
  add_subdirectory(matviz/hom_rot_cross)
  add_subdirectory(matviz/lattice)
  # tests python-matlab interface via doctest
  #add_subdirectory(matlab-interface)
  
  ### commented test as they consume much time
  #add_subdirectory(matviz/interpretation_2x2x2)
  #add_subdirectory(matviz/interpretation_2x2x2_smoothed)
  #add_subdirectory(basecell/interpolation_linear/small_stiffnesses)
  #add_subdirectory(basecell/interpolation_linear/equal_yzstiffnesses)
  #add_subdirectory(basecell/interpolation_linear/equal_stiffnesses)
  #add_subdirectory(basecell/interpolation_linear/surface_mesh/equal_xystiffnesses)
  #add_subdirectory(basecell/interpolation_linear/surface_mesh/equal_stiffnesses_05)
  #add_subdirectory(basecell/interpolation_linear/surface_mesh/equal_stiffnesses_01_res_40)
  #add_subdirectory(basecell/interpolation_linear/surface_mesh/equal_stiffnesses_01_res_20)
  #add_subdirectory(basecell/interpolation_heaviside/small_stiffnesses)
  #add_subdirectory(basecell/interpolation_heaviside/equal_yzstiffnesses)
  #add_subdirectory(basecell/interpolation_heaviside/equal_stiffnesses)
  #add_subdirectory(basecell/interpolation_heaviside/surface_mesh/equal_xystiffnesses)
  #add_subdirectory(basecell/interpolation_heaviside/surface_mesh/equal_stiffnesses_05)
  #add_subdirectory(basecell/interpolation_heaviside/surface_mesh/equal_stiffnesses_01_res_40)
  #add_subdirectory(basecell/interpolation_heaviside/surface_mesh/equal_stiffnesses_01_res_20)
endif(TESTSUITE_PYTHON) 
