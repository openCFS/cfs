# Python binding for muParser

option(BUILD_PYMUPARSER "Build Python binding for muParser" OFF)
mark_as_advanced(BUILD_PYMUPARSER)

if(BUILD_PYMUPARSER AND NOT Python_EXECUTABLE)
    message(WARNING "Cannot build PyMuParser without a Python interpreter.\nPlease set Python_EXECUTABLE")
endif()

if(BUILD_PYMUPARSER AND Python_EXECUTABLE)

    set(PYMUPARSER_VERSION "0.1.3")
    
    configure_file( setup.py.in
                    "${CMAKE_BINARY_DIR}/PyMuParser/setup.py"
                    @ONLY )
    configure_file( muparsermodule.cpp.in
                    "${CMAKE_BINARY_DIR}/PyMuParser/muparsermodule.cpp"
                    @ONLY ) 
    
    add_custom_target(pymuparser ALL
                      ${Python_EXECUTABLE} setup.py build_ext
                      COMMAND ${Python_EXECUTABLE} setup.py bdist_wheel
                      DEPENDS muparser
                      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/PyMuParser"
                      COMMENT "Building PyMuParser wheel")
    
    add_custom_target(install-pymuparser
                      ${Python_EXECUTABLE} -m pip install "dist/PyMuParser-${PYMUPARSER_VERSION}-*.whl"
                      DEPENDS pymuparser
                      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/PyMuParser"
                      COMMENT "Installing PyMuParser wheel")

endif(BUILD_PYMUPARSER AND Python_EXECUTABLE)
