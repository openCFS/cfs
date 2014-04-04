add_external_project(
  ffmpeg
  CONFIGURE_COMMAND <SOURCE_DIR>/configure
                    --prefix=<INSTALL_DIR>
                    --disable-avdevice
                    --disable-bzlib
                    --disable-decoders
                    --disable-doc
                    --disable-ffplay
                    --disable-ffprobe
                    --disable-ffserver
                    --disable-network
                    --disable-static
                    --enable-shared
#                    --extra-cflags=${cppflags}
#                    --extra-ldflags=${ldflags}
  BUILD_IN_SOURCE 1
)
