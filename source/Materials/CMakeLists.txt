SET(MATERIALS_SRCS
  AcousticMaterial.cc
  BaseMaterial.cc
  ElectroMagneticMaterial.cc 
  ElectroStaticMaterial.cc 
  ElecQuasistaticMaterial.cc
  HeatMaterial.cc
  MechanicMaterial.cc
  SmoothMaterial.cc
  FlowMaterial.cc 
  PiezoMaterial.cc
  MagStrictMaterial.cc
  ElectricConductionMaterial.cc
  TestMaterial.cc
  Models/Hysteresis.cc
  Models/EBHysteresis.cc
  Models/SMSM.cc
  Models/Jiles.cc
  Models/Model.cc
  Models/PiezoMicroModelBK.cc
  Models/PiezoMicroModelHF.cc
  Models/VectorPreisachSutor.cc
  Models/Preisach.cc
  Models/VectorPreisachMayergoyz.cc
  Models/SimplePreisachInv.cc 
  )



ADD_LIBRARY(materials STATIC ${MATERIALS_SRCS})

set_target_properties(materials PROPERTIES CXX_STANDARD 17)
set_target_properties(materials PROPERTIES CXX_STANDARD_REQUIRED ON)

# Add necessary include directories for LibTorch
IF(USE_LIBTORCH)
  message(STATUS "LibTorch is enabled and will be included in the build.")
  message(STATUS "LibTorch include directory: ${LIBTORCH_INCLUDE_DIR}")
  message(STATUS "LibTorch libraries: ${LIBTORCH_LIBRARY}")
  INCLUDE_DIRECTORIES(${LIBTORCH_INCLUDE_DIR})  # Ensure LIBTORCH_INCLUDE_DIRS is properly set in your global configuration
  set(LDFLAGS "-L/opt/homebrew/opt/libomp/lib")
  SET(LIBTORCH_LIBRARIES 
  # "/Users/kroppert/Devel/CFS_BIN/buildopt/cfsdeps/libtorch/src/libtorch/lib/*.dylib"
  "/Users/kroppert/Devel/CFS_BIN/buildopt/cfsdeps/libtorch/src/libtorch/lib/libtorch.dylib"
  "/Users/kroppert/Devel/CFS_BIN/buildopt/cfsdeps/libtorch/src/libtorch/lib/libc10.dylib"
  "/Users/kroppert/Devel/CFS_BIN/buildopt/cfsdeps/libtorch/src/libtorch/lib/libtorch_cpu.dylib"
  )


  set(LDFLAGS "-L/opt/homebrew/opt/libomp/lib")
  # Print each file to see if the files are correctly being detected
  file(GLOB TORCH_DYLIBS "/Users/kroppert/Devel/CFS_BIN/buildopt/cfsdeps/libtorch/src/libtorch/lib/*.dylib")
  file(GLOB LIBTORCH_CONTENTS "/Users/kroppert/Devel/CFS_BIN/buildopt/cfsdeps/libtorch/src/libtorch/lib/*")
  message(STATUS "LibTorch directory contents: ${LIBTORCH_CONTENTS}")

  set(MAT_LIB_DIR "$<TARGET_FILE_DIR:materials>")  # Destination directory for copying .... worked

  message(STATUS "DYLIB files found: ${TORCH_DYLIBS}")

  foreach(dylib ${TORCH_DYLIBS})
    message(STATUS "Copying ${dylib} to ${MAT_LIB_DIR}")
    add_custom_command(TARGET materials
                        POST_BUILD
                        COMMAND /bin/bash -c "cp ${dylib} ${MAT_LIB_DIR}/"
                        COMMENT "Copying ${dylib} to ${MAT_LIB_DIR}")
  endforeach()

  # file(GLOB TORCH_DYLIBS "$/Users/kroppert/Devel/CFS_BIN/buildopt/cfsdeps/libtorch/src/libtorch/lib/*.dylib")
  # message("DYLIB files: ${TORCH_DYLIBS}")
  # set(MAT_LIB_DIR "$<TARGET_FILE_DIR:materials>")  # Destination directory for copying
  # message("MAT_LIB_DIR : ${MAT_LIB_DIR}")
  # add_custom_command(TARGET materials
  #                    POST_BUILD
  #                    COMMAND cp ${TORCH_DYLIBS} ${MAT_LIB_DIR}
  #                    COMMENT "Copying LibTorch dylibs to ${MAT_LIB_DIR}")
  
  
#   SET(TARGET_LL ${TARGET_LL} libtorch-materials)

ENDIF(USE_LIBTORCH)

TARGET_LINK_LIBRARIES(materials
  utils
  odesolver
  coeffunction
  ${LIBTORCH_LIBRARIES}
)

