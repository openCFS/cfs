# default Windows setting overwriting platform_defaults.cmake 
# will be overrided e.g. by a .cfs_platform_defaults.cmake in the user's home dire
set(USE_GIDPOST_DEFAULT OFF)
set(USE_SUPERLU_DEFAULT OFF)
set(USE_CGNS_DEFAULT OFF)
set(USE_GINKGO_DEFAULT ON) # requires clang (e.g. from MSVC)

