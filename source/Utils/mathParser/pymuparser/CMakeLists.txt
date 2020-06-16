# Python binding for muParser

option(BUILD_PYMUPARSER "Build and install Python binding for muParser" OFF)
mark_as_advanced(BUILD_PYMUPARSER)

if(NOT PYTHON_EXECUTABLE)
    message(WARNING "Cannot build PyMuParser without a Python interpreter.\nPlease set PYTHON_EXECUTABLE")
endif(NOT PYTHON_EXECUTABLE)

if(BUILD_PYMUPARSER AND PYTHON_EXECUTABLE)

    set(PYMUPARSER_VERSION "0.1.2")
    
    configure_file( setup.py.in
                    "${CMAKE_BINARY_DIR}/PyMuParser/setup.py"
                    @ONLY )
    configure_file( muparsermodule.cpp.in
                    "${CMAKE_BINARY_DIR}/PyMuParser/muparsermodule.cpp"
                    @ONLY ) 
    
    add_custom_target(PyMuParser ALL
                      ${PYTHON_EXECUTABLE} setup.py build_ext
                      COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_wheel
                      DEPENDS muparser
                      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/PyMuParser"
                      COMMENT "Building PyMuParser wheel")
    
endif(BUILD_PYMUPARSER AND PYTHON_EXECUTABLE)
