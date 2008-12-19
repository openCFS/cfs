SET(SIMOUTVTK_SRCS SimOutputVTK.cc)

# where the original stuff resides
SET(VTK_DIR /home/local/lse/fwein/packages/ParaView-3.4.0)

# where we build
SET(VTK_BUILD ${VTK_DIR}/build)

# abbreviation
SET(VTK_INCLUDE_BASE ${VTK_DIR}/VTK)
INCLUDE_DIRECTORIES(
    ${VTK_BUILD}/VTK
    ${VTK_INCLUDE_BASE}/IO
    ${VTK_INCLUDE_BASE}/Filtering
    ${VTK_INCLUDE_BASE}/Common
  )
    
ADD_LIBRARY(simoutvtk STATIC ${SIMOUTVTK_SRCS})

# add this also to LD_LIBRARY_PATH!
TARGET_LINK_LIBRARIES(simoutvtk
    ${VTK_DIR}/build/bin/libvtkIO.so
    ${VTK_DIR}/build/bin/libvtkFiltering.so
    ${VTK_DIR}/build/bin/libvtkCommon.so
)


