add_external_project(
  freetype
  CONFIGURE_COMMAND <SOURCE_DIR>/configure
                    --prefix=<INSTALL_DIR>
                    --enable-static=no
                    --with-sysroot=<INSTALL_DIR>
) 
